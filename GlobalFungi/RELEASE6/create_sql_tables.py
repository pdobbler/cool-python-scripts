#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from contextlib import ExitStack


MAX_UINT32 = 4294967295


def parse_arguments():

    parser = argparse.ArgumentParser(
        description=(
            "Convert a variant TSV table into four MariaDB "
            "import tables: variants, samplevar, sh_table "
            "and sample_table."
        )
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input variants TSV file"
    )

    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files"
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write column headers"
    )

    return parser.parse_args()


def prepare_output_files(outdir, overwrite):

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = {
        "variants": outdir / "variants.tsv",
        "samplevar": outdir / "samplevar.tsv",
        "sh_table": outdir / "sh_table.tsv",
        "sample_table": outdir / "sample_table.tsv"
    }

    if not overwrite:

        existing = [
            str(path)
            for path in paths.values()
            if path.exists()
        ]

        if existing:
            raise FileExistsError(
                "Output files already exist:\n"
                + "\n".join(existing)
                + "\nUse --overwrite to replace them."
            )

    return paths


def validate_header(header):

    required = [
        "seqID",
        "samples",
        "abundances",
        "HIT",
        "similarity",
        "coverage",
        "marker",
        "sequence"
    ]

    if len(header) != len(set(header)):
        raise ValueError(
            "Input contains duplicate column names."
        )

    missing = [
        column
        for column in required
        if column not in header
    ]

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
        )

    return {
        name: header.index(name)
        for name in required
    }


def main():

    args = parse_arguments()

    input_file = Path(args.input)

    if not input_file.is_file():
        raise FileNotFoundError(
            f"Input file does not exist: {input_file}"
        )

    paths = prepare_output_files(
        args.outdir,
        args.overwrite
    )

    print(
        f"Reading input: {input_file}",
        file=sys.stderr
    )

    # ---------------------------------------------------------
    # ID dictionaries stored in RAM
    # ---------------------------------------------------------

    sample_ids = {}
    sh_ids = {}

    # ---------------------------------------------------------
    # Counters
    # ---------------------------------------------------------

    variant_id = 0
    samplevar_id = 0

    no_hit_variants = 0
    identified_variants = 0

    total_abundance = 0

    # ---------------------------------------------------------
    # Open input and all output files
    # ---------------------------------------------------------

    with ExitStack() as stack:

        infile = stack.enter_context(
            input_file.open(
                "r",
                encoding="utf-8",
                buffering=1024 * 1024
            )
        )

        variants_out = stack.enter_context(
            paths["variants"].open(
                "w",
                encoding="utf-8",
                buffering=1024 * 1024
            )
        )

        samplevar_out = stack.enter_context(
            paths["samplevar"].open(
                "w",
                encoding="utf-8",
                buffering=1024 * 1024
            )
        )

        sh_out = stack.enter_context(
            paths["sh_table"].open(
                "w",
                encoding="utf-8",
                buffering=1024 * 1024
            )
        )

        sample_out = stack.enter_context(
            paths["sample_table"].open(
                "w",
                encoding="utf-8",
                buffering=1024 * 1024
            )
        )

        # -----------------------------------------------------
        # Read input header
        # -----------------------------------------------------

        header_line = infile.readline()

        if not header_line:
            raise ValueError("Input file is empty.")

        header = header_line.rstrip("\r\n").split("\t")

        col = validate_header(header)

        expected_fields = len(header)

        # -----------------------------------------------------
        # Write output headers
        # -----------------------------------------------------

        if not args.no_header:

            variants_out.write(
                "id\tSH_id\tmarker\thash\tsequence\n"
            )

            samplevar_out.write(
                "id\tvariant\tsample\tabundance\tsh_id\n"
            )

            sh_out.write(
                "SH_id\tSH_name\n"
            )

            sample_out.write(
                "sample\tsample_name\n"
            )

        # -----------------------------------------------------
        # Process input
        # -----------------------------------------------------

        for line_number, line in enumerate(
            infile,
            start=2
        ):

            if not line.strip():
                continue

            fields = line.rstrip("\r\n").split("\t")

            if len(fields) != expected_fields:

                raise ValueError(
                    f"Line {line_number}: expected "
                    f"{expected_fields} columns, "
                    f"found {len(fields)}."
                )

            seq_id = fields[col["seqID"]]
            samples_text = fields[col["samples"]]
            abundances_text = fields[col["abundances"]]
            hit = fields[col["HIT"]].strip()
            marker = fields[col["marker"]]
            sequence = fields[col["sequence"]]

            # -------------------------------------------------
            # Basic validation
            # -------------------------------------------------

            if len(seq_id) != 32:

                raise ValueError(
                    f"Line {line_number}: invalid seqID: "
                    f"{seq_id}"
                )

            if not marker or len(marker) > 4:

                raise ValueError(
                    f"Line {line_number}: marker must "
                    f"contain 1-4 characters: {marker}"
                )

            if not sequence:

                raise ValueError(
                    f"Line {line_number}: empty sequence."
                )

            if not samples_text or not abundances_text:

                raise ValueError(
                    f"Line {line_number}: missing samples "
                    f"or abundances."
                )

            if not hit:

                raise ValueError(
                    f"Line {line_number}: empty HIT."
                )

            # -------------------------------------------------
            # Parse samples and abundances
            # -------------------------------------------------

            samples = samples_text.split(";")
            abundances = abundances_text.split(";")

            if len(samples) != len(abundances):

                raise ValueError(
                    f"Line {line_number}: number of samples "
                    f"({len(samples)}) does not match number "
                    f"of abundances ({len(abundances)})."
                )

            if len(samples) != len(set(samples)):

                raise ValueError(
                    f"Line {line_number}: duplicated sample "
                    f"within a single variant."
                )

            parsed_abundances = []

            for abundance_text in abundances:

                try:
                    abundance = int(abundance_text)

                except ValueError:
                    raise ValueError(
                        f"Line {line_number}: invalid abundance: "
                        f"{abundance_text}"
                    )

                if abundance < 1 or abundance > MAX_UINT32:

                    raise ValueError(
                        f"Line {line_number}: abundance outside "
                        f"MariaDB unsigned INT range: {abundance}"
                    )

                parsed_abundances.append(abundance)

            # -------------------------------------------------
            # Assign variant ID
            # -------------------------------------------------

            variant_id += 1

            if variant_id > MAX_UINT32:

                raise OverflowError(
                    "Variant ID exceeds MariaDB unsigned INT."
                )

            # -------------------------------------------------
            # Assign SH_id
            # -------------------------------------------------

            # No identification:
            # HIT = "-" or "NA" -> SH_id = 0

            if hit in ("-", "NA"):

                sh_id = 0
                no_hit_variants += 1

            else:

                identified_variants += 1

                if hit not in sh_ids:

                    new_sh_id = len(sh_ids) + 1

                    if new_sh_id > MAX_UINT32:

                        raise OverflowError(
                            "SH_id exceeds MariaDB unsigned INT."
                        )

                    sh_ids[hit] = new_sh_id

                    # Write SH mapping immediately
                    sh_out.write(
                        f"{new_sh_id}\t{hit}\n"
                    )

                sh_id = sh_ids[hit]

            # -------------------------------------------------
            # Write variants table
            # -------------------------------------------------

            variants_out.write(
                f"{variant_id}\t"
                f"{sh_id}\t"
                f"{marker}\t"
                f"{seq_id}\t"
                f"{sequence}\n"
            )

            # -------------------------------------------------
            # Process individual samples
            # -------------------------------------------------

            for sample_name, abundance in zip(
                samples,
                parsed_abundances
            ):

                if not sample_name:

                    raise ValueError(
                        f"Line {line_number}: empty sample name."
                    )

                # ---------------------------------------------
                # Assign sample ID
                # ---------------------------------------------

                if sample_name not in sample_ids:

                    new_sample_id = len(sample_ids) + 1

                    if new_sample_id > MAX_UINT32:

                        raise OverflowError(
                            "Sample ID exceeds MariaDB unsigned INT."
                        )

                    sample_ids[sample_name] = new_sample_id

                    # Write mapping immediately
                    sample_out.write(
                        f"{new_sample_id}\t"
                        f"{sample_name}\n"
                    )

                sample_id = sample_ids[sample_name]

                # ---------------------------------------------
                # Assign samplevar ID
                # ---------------------------------------------

                samplevar_id += 1

                # samplevar.id is unsigned BIGINT.
                # Python integers support this range.

                if samplevar_id > 18446744073709551615:

                    raise OverflowError(
                        "samplevar ID exceeds unsigned BIGINT."
                    )

                # ---------------------------------------------
                # Write samplevar record
                # ---------------------------------------------

                samplevar_out.write(
                    f"{samplevar_id}\t"
                    f"{variant_id}\t"
                    f"{sample_id}\t"
                    f"{abundance}\t"
                    f"{sh_id}\n"
                )

                total_abundance += abundance

            # -------------------------------------------------
            # Progress report every 1 million variants
            # -------------------------------------------------

            if variant_id % 1_000_000 == 0:

                print(
                    f"Processed variants: {variant_id:,} | "
                    f"samplevar records: {samplevar_id:,} | "
                    f"samples: {len(sample_ids):,} | "
                    f"SH IDs: {len(sh_ids):,}",
                    file=sys.stderr
                )

    # ---------------------------------------------------------
    # Final report
    # ---------------------------------------------------------

    print(
        "\n=== FINAL REPORT ===",
        file=sys.stderr
    )

    print(
        f"Variants:                 {variant_id:,}",
        file=sys.stderr
    )

    print(
        f"Samplevar records:        {samplevar_id:,}",
        file=sys.stderr
    )

    print(
        f"Unique samples:           {len(sample_ids):,}",
        file=sys.stderr
    )

    print(
        f"Unique SH IDs:            {len(sh_ids):,}",
        file=sys.stderr
    )

    print(
        f"Identified variants:      {identified_variants:,}",
        file=sys.stderr
    )

    print(
        f"Variants without HIT:     {no_hit_variants:,}",
        file=sys.stderr
    )

    print(
        f"Total abundance:          {total_abundance:,}",
        file=sys.stderr
    )

    print(
        "\nOutput files:",
        file=sys.stderr
    )

    for name, path in paths.items():

        print(
            f"  {name}: {path}",
            file=sys.stderr
        )

    print(
        "\nFinished successfully.",
        file=sys.stderr
    )


if __name__ == "__main__":

    try:
        main()

    except (OSError, ValueError, OverflowError) as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr
        )

        sys.exit(1)
