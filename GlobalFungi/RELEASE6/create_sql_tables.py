#!/usr/bin/env python3

import argparse
import gzip
import re
import sys

from pathlib import Path
from contextlib import ExitStack


MAX_UINT32 = 4294967295
MAX_UINT64 = 18446744073709551615

SEQID_PATTERN = re.compile(r"[0-9a-fA-F]{32}")


# =============================================================
# ARGUMENTS
# =============================================================

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
        help="Input variants TSV file (plain text or .gz)"
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
        help="Do not write output column headers"
    )

    return parser.parse_args()


# =============================================================
# INPUT
# =============================================================

def open_input(filename):

    if str(filename).endswith(".gz"):

        return gzip.open(
            filename,
            "rt",
            encoding="utf-8"
        )

    return open(
        filename,
        "r",
        encoding="utf-8",
        buffering=1024 * 1024
    )


# =============================================================
# OUTPUT FILES
# =============================================================

def prepare_output_files(outdir, overwrite, input_file):

    outdir = Path(outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    paths = {
        "variants": outdir / "variants.tsv",
        "samplevar": outdir / "samplevar.tsv",
        "sh_table": outdir / "sh_table.tsv",
        "sample_table": outdir / "sample_table.tsv"
    }

    # Prevent accidental overwriting of input file.

    for path in paths.values():

        if path.resolve() == input_file.resolve():

            raise ValueError(
                f"Input file cannot also be an output file: {path}"
            )

    # Check whether output files already exist.

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


# =============================================================
# VALIDATE HEADER
# =============================================================

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
            "Input contains duplicated column names."
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


# =============================================================
# MAIN
# =============================================================

def main():

    args = parse_arguments()

    input_file = Path(args.input)

    if not input_file.is_file():

        raise FileNotFoundError(
            f"Input file does not exist: {input_file}"
        )

    paths = prepare_output_files(
        args.outdir,
        args.overwrite,
        input_file
    )

    print(
        f"Reading input: {input_file}",
        file=sys.stderr
    )

    # =========================================================
    # DICTIONARIES IN RAM
    # =========================================================

    # sample_name -> numerical sample ID

    sample_ids = {}

    # SH_name -> numerical SH ID

    sh_ids = {}

    # Previously encountered sequence IDs.
    # Store 16-byte MD5 digests rather than 32-character strings
    # to reduce memory usage.

    seen_seqids = set()

    # =========================================================
    # DUPLICATE seqID STATISTICS
    # =========================================================

    duplicate_seqid_count = 0

    duplicate_seqid_examples = []

    MAX_DUPLICATE_EXAMPLES = 20

    # =========================================================
    # COUNTERS
    # =========================================================

    variant_id = 0

    samplevar_id = 0

    identified_variants = 0

    no_hit_variants = 0

    total_abundance = 0

    # =========================================================
    # OPEN FILES
    # =========================================================

    with ExitStack() as stack:

        infile = stack.enter_context(
            open_input(input_file)
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

        # =====================================================
        # READ INPUT HEADER
        # =====================================================

        header_line = infile.readline()

        if not header_line:

            raise ValueError(
                "Input file is empty."
            )

        header = header_line.rstrip("\r\n").split("\t")

        col = validate_header(header)

        expected_fields = len(header)

        # =====================================================
        # WRITE OUTPUT HEADERS
        # =====================================================

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

        # =====================================================
        # PROCESS INPUT TABLE
        # =====================================================

        for line_number, line in enumerate(
            infile,
            start=2
        ):

            if not line.strip():
                continue

            fields = line.rstrip("\r\n").split("\t")

            # -------------------------------------------------
            # VALIDATE COLUMN COUNT
            # -------------------------------------------------

            if len(fields) != expected_fields:

                raise ValueError(
                    f"Line {line_number}: expected "
                    f"{expected_fields} columns, "
                    f"found {len(fields)}."
                )

            # -------------------------------------------------
            # EXTRACT VALUES
            # -------------------------------------------------

            seq_id = fields[col["seqID"]]

            samples_text = fields[col["samples"]]

            abundances_text = fields[col["abundances"]]

            hit = fields[col["HIT"]].strip()

            marker = fields[col["marker"]]

            sequence = fields[col["sequence"]]

            # =================================================
            # BASIC VALIDATION
            # =================================================

            # -------------------------------------------------
            # Validate seqID
            # -------------------------------------------------

            if not SEQID_PATTERN.fullmatch(seq_id):

                raise ValueError(
                    f"Line {line_number}: invalid seqID: "
                    f"{seq_id}"
                )

            # -------------------------------------------------
            # Validate marker
            # -------------------------------------------------

            if not marker or len(marker) > 4:

                raise ValueError(
                    f"Line {line_number}: marker must "
                    f"contain 1-4 characters: {marker}"
                )

            # -------------------------------------------------
            # Validate sequence
            # -------------------------------------------------

            if not sequence:

                raise ValueError(
                    f"Line {line_number}: empty sequence."
                )

            # -------------------------------------------------
            # Validate samples and abundances
            # -------------------------------------------------

            if not samples_text or not abundances_text:

                raise ValueError(
                    f"Line {line_number}: missing samples "
                    f"or abundances."
                )

            # -------------------------------------------------
            # Validate HIT
            # -------------------------------------------------

            if not hit:

                raise ValueError(
                    f"Line {line_number}: empty HIT."
                )

            # =================================================
            # PARSE SAMPLES AND ABUNDANCES
            # =================================================

            samples = samples_text.split(";")

            abundances = abundances_text.split(";")

            # -------------------------------------------------
            # Check matching sample/abundance counts
            # -------------------------------------------------

            if len(samples) != len(abundances):

                raise ValueError(
                    f"Line {line_number}: number of samples "
                    f"({len(samples)}) does not match number "
                    f"of abundances ({len(abundances)})."
                )

            # -------------------------------------------------
            # Check for duplicated samples within variant
            # -------------------------------------------------

            if len(samples) != len(set(samples)):

                raise ValueError(
                    f"Line {line_number}: duplicated sample "
                    f"within a single variant."
                )

            # -------------------------------------------------
            # Check empty sample names
            # -------------------------------------------------

            for sample_name in samples:

                if not sample_name:

                    raise ValueError(
                        f"Line {line_number}: empty sample name."
                    )

            # -------------------------------------------------
            # Parse abundances
            # -------------------------------------------------

            parsed_abundances = []

            for abundance_text in abundances:

                try:

                    abundance = int(abundance_text)

                except ValueError:

                    raise ValueError(
                        f"Line {line_number}: invalid abundance: "
                        f"{abundance_text}"
                    )

                # MariaDB unsigned INT range

                if abundance < 1 or abundance > MAX_UINT32:

                    raise ValueError(
                        f"Line {line_number}: abundance outside "
                        f"MariaDB unsigned INT range: {abundance}"
                    )

                parsed_abundances.append(abundance)

            # =================================================
            # CHECK DUPLICATE seqID
            # =================================================

            # Convert MD5 hex string into 16-byte representation.
            # This reduces RAM consumption compared with storing
            # the original 32-character strings.

            seq_id_binary = bytes.fromhex(seq_id)

            if seq_id_binary in seen_seqids:

                duplicate_seqid_count += 1

                # Keep only first 20 examples.

                if (
                    len(duplicate_seqid_examples)
                    < MAX_DUPLICATE_EXAMPLES
                ):

                    duplicate_seqid_examples.append(
                        (seq_id, line_number)
                    )

            else:

                seen_seqids.add(seq_id_binary)

            # IMPORTANT:
            # Duplicated seqIDs are reported but NOT skipped.

            # =================================================
            # ASSIGN VARIANT ID
            # =================================================

            variant_id += 1

            if variant_id > MAX_UINT32:

                raise OverflowError(
                    "Variant ID exceeds MariaDB unsigned INT."
                )

            # =================================================
            # ASSIGN SH_id
            # =================================================

            # HIT = "-" or "NA" -> SH_id = 0
            # Both variants and samplevar use the same SH_id.

            if hit in ("-", "NA"):

                sh_id = 0

                no_hit_variants += 1

            else:

                identified_variants += 1

                # ---------------------------------------------
                # New SH name
                # ---------------------------------------------

                if hit not in sh_ids:

                    new_sh_id = len(sh_ids) + 1

                    if new_sh_id > MAX_UINT32:

                        raise OverflowError(
                            "SH_id exceeds MariaDB unsigned INT."
                        )

                    sh_ids[hit] = new_sh_id

                    # Write SH mapping.

                    sh_out.write(
                        f"{new_sh_id}\t{hit}\n"
                    )

                # Get assigned SH ID.

                sh_id = sh_ids[hit]

            # =================================================
            # WRITE VARIANTS TABLE
            # =================================================

            variants_out.write(
                f"{variant_id}\t"
                f"{sh_id}\t"
                f"{marker}\t"
                f"{seq_id}\t"
                f"{sequence}\n"
            )

            # =================================================
            # PROCESS INDIVIDUAL SAMPLES
            # =================================================

            for sample_name, abundance in zip(
                samples,
                parsed_abundances
            ):

                # =============================================
                # ASSIGN SAMPLE ID
                # =============================================

                if sample_name not in sample_ids:

                    new_sample_id = len(sample_ids) + 1

                    if new_sample_id > MAX_UINT32:

                        raise OverflowError(
                            "Sample ID exceeds MariaDB unsigned INT."
                        )

                    sample_ids[sample_name] = new_sample_id

                    # Write sample mapping.

                    sample_out.write(
                        f"{new_sample_id}\t"
                        f"{sample_name}\n"
                    )

                sample_id = sample_ids[sample_name]

                # =============================================
                # ASSIGN SAMPLEVAR ID
                # =============================================

                samplevar_id += 1

                if samplevar_id > MAX_UINT64:

                    raise OverflowError(
                        "samplevar ID exceeds unsigned BIGINT."
                    )

                # =============================================
                # WRITE SAMPLEVAR RECORD
                # =============================================

                samplevar_out.write(
                    f"{samplevar_id}\t"
                    f"{variant_id}\t"
                    f"{sample_id}\t"
                    f"{abundance}\t"
                    f"{sh_id}\n"
                )

                # =============================================
                # UPDATE TOTAL ABUNDANCE
                # =============================================

                total_abundance += abundance

            # =================================================
            # PROGRESS REPORT
            # =================================================

            # Report every 1 million processed variants.

            if variant_id % 1_000_000 == 0:

                print(
                    f"Processed variants: {variant_id:,} | "
                    f"samplevar records: {samplevar_id:,} | "
                    f"samples: {len(sample_ids):,} | "
                    f"SH IDs: {len(sh_ids):,} | "
                    f"duplicate seqIDs: {duplicate_seqid_count:,}",
                    file=sys.stderr
                )

    # =========================================================
    # FINAL REPORT
    # =========================================================

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

    # =========================================================
    # DUPLICATE seqID REPORT
    # =========================================================

    print(
        f"Duplicate seqID occurrences: {duplicate_seqid_count:,}",
        file=sys.stderr
    )

    print(
        f"Unique seqIDs:            {len(seen_seqids):,}",
        file=sys.stderr
    )

    if duplicate_seqid_count > 0:

        print(
            "\nWARNING: Duplicate seqIDs detected!",
            file=sys.stderr
        )

        print(
            "Duplicate records were NOT skipped.",
            file=sys.stderr
        )

        print(
            f"First {MAX_DUPLICATE_EXAMPLES} "
            f"duplicate occurrences:",
            file=sys.stderr
        )

        for seq_id, line_number in duplicate_seqid_examples:

            print(
                f"  Line {line_number}: {seq_id}",
                file=sys.stderr
            )

        if duplicate_seqid_count > MAX_DUPLICATE_EXAMPLES:

            print(
                f"  ... and "
                f"{duplicate_seqid_count - MAX_DUPLICATE_EXAMPLES:,} "
                f"more duplicate occurrences",
                file=sys.stderr
            )

    else:

        print(
            "\nNo duplicate seqIDs detected.",
            file=sys.stderr
        )

    # =========================================================
    # OUTPUT FILE REPORT
    # =========================================================

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


# =============================================================
# EXECUTION
# =============================================================

if __name__ == "__main__":

    try:

        main()

    except (
        OSError,
        ValueError,
        OverflowError
    ) as error:

        print(
            f"\nERROR: {error}",
            file=sys.stderr
        )

        print(
            "Processing failed. Output files may be incomplete.",
            file=sys.stderr
        )

        sys.exit(1)
