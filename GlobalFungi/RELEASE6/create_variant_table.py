#!/usr/bin/env python3

import argparse
import gzip
import hashlib
import sys
from collections import Counter, OrderedDict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Create a tab-delimited sequence variant table from a gzipped FASTA "
            "file and a gzipped sequence identification table."
        )
    )

    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Input gzipped FASTA file"
    )

    parser.add_argument(
        "-t", "--table",
        required=True,
        help="Input gzipped identification table"
    )

    parser.add_argument(
        "-s", "--similarity",
        required=True,
        type=float,
        help="Minimum similarity threshold"
    )

    parser.add_argument(
        "-c", "--coverage",
        required=True,
        type=float,
        help="Minimum coverage threshold"
    )

    parser.add_argument(
        "-m", "--marker",
        required=True,
        help="Marker name, e.g. ITS1"
    )

    parser.add_argument(
        "-o", "--output",
        default="-",
        help="Output TSV file. Default: stdout"
    )

    return parser.parse_args()


def open_text_gz(filename):
    """Open gzip-compressed text file."""
    return gzip.open(filename, "rt", encoding="utf-8")


def read_fasta(filename):
    """
    Read gzipped FASTA and aggregate identical sequences.

    For each sequence:
      - sample name is taken from FASTA header before the first "|"
      - seqID is calculated as md5(sequence)
      - abundance is counted separately for each sample

    Returns
    -------
    variants : OrderedDict
        seqID -> {
            "sequence": sequence,
            "samples": Counter({
                sample1: count,
                sample2: count,
                ...
            })
        }

    fasta_records : int
        Total number of FASTA records.
    """

    variants = OrderedDict()

    header = None
    sequence_parts = []
    fasta_records = 0

    def process_record(record_header, seq_parts):
        nonlocal fasta_records

        if record_header is None:
            return

        sequence = "".join(seq_parts)

        if not sequence:
            print(
                f"WARNING: Empty sequence for FASTA header: {record_header}",
                file=sys.stderr
            )
            return

        # Sample name = everything before first "|"
        sample = record_header.split("|", 1)[0]

        # seqID calculated exactly from sequence
        seq_id = hashlib.md5(sequence.encode()).hexdigest()

        if seq_id not in variants:
            variants[seq_id] = {
                "sequence": sequence,
                "samples": Counter()
            }

        else:
            # Extremely unlikely, but check for MD5 collision
            if variants[seq_id]["sequence"] != sequence:
                raise RuntimeError(
                    f"MD5 collision detected for seqID {seq_id}"
                )

        variants[seq_id]["samples"][sample] += 1
        fasta_records += 1

    with open_text_gz(filename) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Process previous record
                process_record(header, sequence_parts)

                header = line[1:]
                sequence_parts = []

            else:
                if header is None:
                    raise ValueError(
                        "FASTA sequence encountered before the first header."
                    )

                sequence_parts.append(line)

        # Process final FASTA record
        process_record(header, sequence_parts)

    return variants, fasta_records


def read_identification_table(filename):
    """
    Read gzipped identification table.

    Expected columns include:

        QUERY
        HIT
        SIMILARITY
        COVERAGE

    Special cases:

    1) Normal HIT:
       similarity and coverage are parsed as floats.

    2) No HIT:
       HIT, SIMILARITY or COVERAGE contains "-"
       -> stored internally as no_hit=True

    3) Duplicate QUERY:
       the first occurrence is used and duplicates are reported.

    Returns
    -------
    identifications : dict

    duplicate_queries : Counter

    table_rows : int
        Number of non-empty data rows.

    no_hit_rows : int
        Number of unique QUERY entries without HIT.
    """

    identifications = {}
    query_counts = Counter()

    table_rows = 0
    no_hit_rows = 0

    required_columns = {
        "QUERY",
        "HIT",
        "SIMILARITY",
        "COVERAGE"
    }

    with open_text_gz(filename) as handle:
        header_line = handle.readline()

        if not header_line:
            raise ValueError(
                "Identification table is empty."
            )

        header = header_line.strip().split()

        missing_columns = required_columns - set(header)

        if missing_columns:
            raise ValueError(
                "Missing required column(s) in identification table: "
                + ", ".join(sorted(missing_columns))
            )

        col = {
            name: index
            for index, name in enumerate(header)
        }

        max_required_index = max(
            col["QUERY"],
            col["HIT"],
            col["SIMILARITY"],
            col["COVERAGE"]
        )

        for line_number, line in enumerate(handle, start=2):
            line = line.strip()

            if not line:
                continue

            table_rows += 1

            fields = line.split()

            if len(fields) <= max_required_index:
                print(
                    f"WARNING: Skipping malformed line "
                    f"{line_number} in identification table:",
                    file=sys.stderr
                )
                print(
                    f"  {line}",
                    file=sys.stderr
                )
                continue

            query = fields[col["QUERY"]]
            hit = fields[col["HIT"]]
            similarity_text = fields[col["SIMILARITY"]]
            coverage_text = fields[col["COVERAGE"]]

            query_counts[query] += 1

            # If QUERY occurs more than once,
            # keep only the first occurrence.
            if query in identifications:
                continue

            # ---------------------------------------------------------
            # NO HIT
            # ---------------------------------------------------------
            #
            # Examples:
            #
            # QUERY   HIT   SIMILARITY   COVERAGE
            # xxx     -     -            -
            #
            # If any of these values is "-", treat the complete
            # identification as absent.
            #
            if (
                hit == "-"
                or similarity_text == "-"
                or coverage_text == "-"
            ):
                identifications[query] = {
                    "hit": "-",
                    "similarity": None,
                    "similarity_text": "-",
                    "coverage": None,
                    "coverage_text": "-",
                    "no_hit": True
                }

                no_hit_rows += 1
                continue

            # ---------------------------------------------------------
            # NORMAL HIT
            # ---------------------------------------------------------

            try:
                similarity = float(similarity_text)
                coverage = float(coverage_text)

            except ValueError:
                raise ValueError(
                    f"Invalid SIMILARITY or COVERAGE "
                    f"at line {line_number}: "
                    f"{similarity_text}, {coverage_text}"
                )

            identifications[query] = {
                "hit": hit,
                "similarity": similarity,
                "similarity_text": similarity_text,
                "coverage": coverage,
                "coverage_text": coverage_text,
                "no_hit": False
            }

    duplicate_queries = Counter(
        {
            query: count
            for query, count in query_counts.items()
            if count > 1
        }
    )

    return (
        identifications,
        duplicate_queries,
        table_rows,
        no_hit_rows
    )


def write_output(
    variants,
    identifications,
    min_similarity,
    min_coverage,
    marker,
    output_filename
):
    """
    Write output TSV.

    Rules:

    QUERY absent from identification table:
        HIT        NA
        similarity NA
        coverage   NA

    QUERY exists but has no HIT:
        HIT        -
        similarity -
        coverage   -

    QUERY exists, HIT exists, but threshold is not satisfied:
        HIT        -
        similarity -
        coverage   -

    QUERY exists and passes thresholds:
        HIT        original HIT
        similarity original similarity
        coverage   original coverage
    """

    if output_filename == "-":
        out = sys.stdout
        close_output = False

    else:
        out = open(
            output_filename,
            "w",
            encoding="utf-8"
        )
        close_output = True

    missing_queries = []

    passed_threshold = 0
    failed_threshold = 0
    no_hit_count = 0

    try:
        # Output header
        print(
            "\t".join([
                "seqID",
                "samples",
                "abundances",
                "HIT",
                "similarity",
                "coverage",
                "marker",
                "sequence"
            ]),
            file=out
        )

        for seq_id, data in variants.items():

            # Keep samples and abundance in exactly the same order.
            sample_names = list(
                data["samples"].keys()
            )

            samples = ";".join(sample_names)

            abundances = ";".join(
                str(data["samples"][sample])
                for sample in sample_names
            )

            # ---------------------------------------------------------
            # QUERY NOT PRESENT AT ALL
            # ---------------------------------------------------------

            if seq_id not in identifications:

                hit_out = "NA"
                similarity_out = "NA"
                coverage_out = "NA"

                missing_queries.append(seq_id)

            else:

                identification = identifications[seq_id]

                # -----------------------------------------------------
                # QUERY EXISTS, BUT THERE IS NO HIT
                # -----------------------------------------------------

                if identification["no_hit"]:

                    hit_out = "-"
                    similarity_out = "-"
                    coverage_out = "-"

                    no_hit_count += 1

                # -----------------------------------------------------
                # HIT EXISTS AND PASSES THRESHOLDS
                # -----------------------------------------------------

                elif (
                    identification["similarity"] >= min_similarity
                    and
                    identification["coverage"] >= min_coverage
                ):

                    hit_out = identification["hit"]

                    similarity_out = (
                        identification["similarity_text"]
                    )

                    coverage_out = (
                        identification["coverage_text"]
                    )

                    passed_threshold += 1

                # -----------------------------------------------------
                # HIT EXISTS BUT FAILS THRESHOLDS
                # -----------------------------------------------------

                else:

                    hit_out = "-"
                    similarity_out = "-"
                    coverage_out = "-"

                    failed_threshold += 1

            print(
                "\t".join([
                    seq_id,
                    samples,
                    abundances,
                    hit_out,
                    similarity_out,
                    coverage_out,
                    marker,
                    data["sequence"]
                ]),
                file=out
            )

    finally:
        if close_output:
            out.close()

    return {
        "missing_queries": missing_queries,
        "passed_threshold": passed_threshold,
        "failed_threshold": failed_threshold,
        "no_hit_count": no_hit_count
    }


def main():

    args = parse_arguments()

    print(
        f"Reading FASTA: {args.fasta}",
        file=sys.stderr
    )

    variants, fasta_records = read_fasta(
        args.fasta
    )

    print(
        f"FASTA records: {fasta_records}",
        file=sys.stderr
    )

    print(
        f"Unique sequence variants: {len(variants)}",
        file=sys.stderr
    )

    print(
        f"Reading identification table: {args.table}",
        file=sys.stderr
    )

    (
        identifications,
        duplicate_queries,
        table_rows,
        identification_no_hit_rows
    ) = read_identification_table(
        args.table
    )

    print(
        f"Identification table rows: {table_rows}",
        file=sys.stderr
    )

    print(
        f"Unique QUERY IDs in identification table: "
        f"{len(identifications)}",
        file=sys.stderr
    )

    print(
        f"QUERY IDs without HIT: "
        f"{identification_no_hit_rows}",
        file=sys.stderr
    )

    print(
        f"Writing output: {args.output}",
        file=sys.stderr
    )

    statistics = write_output(
        variants=variants,
        identifications=identifications,
        min_similarity=args.similarity,
        min_coverage=args.coverage,
        marker=args.marker,
        output_filename=args.output
    )

    missing_queries = statistics["missing_queries"]

    # =============================================================
    # FINAL REPORT
    # =============================================================

    print(
        "",
        file=sys.stderr
    )

    print(
        "=== FINAL REPORT ===",
        file=sys.stderr
    )

    print(
        f"FASTA records:                  "
        f"{fasta_records}",
        file=sys.stderr
    )

    print(
        f"Unique sequence variants:       "
        f"{len(variants)}",
        file=sys.stderr
    )

    print(
        f"Identification table rows:      "
        f"{table_rows}",
        file=sys.stderr
    )

    print(
        f"Unique identification QUERYs:   "
        f"{len(identifications)}",
        file=sys.stderr
    )

    print(
        f"Passed thresholds:              "
        f"{statistics['passed_threshold']}",
        file=sys.stderr
    )

    print(
        f"Failed thresholds:              "
        f"{statistics['failed_threshold']}",
        file=sys.stderr
    )

    print(
        f"QUERYs without HIT (-):         "
        f"{statistics['no_hit_count']}",
        file=sys.stderr
    )

    print(
        f"FASTA seqIDs without QUERY:     "
        f"{len(missing_queries)}",
        file=sys.stderr
    )

    print(
        f"Duplicated QUERY IDs:           "
        f"{len(duplicate_queries)}",
        file=sys.stderr
    )

    # =============================================================
    # MISSING QUERY WARNING
    # =============================================================

    if missing_queries:

        print(
            "",
            file=sys.stderr
        )

        print(
            "WARNING: Some FASTA seqIDs were not found "
            "as QUERY in the identification table.",
            file=sys.stderr
        )

        print(
            "Their HIT/SIMILARITY/COVERAGE values "
            "were set to NA.",
            file=sys.stderr
        )

        max_show = 20

        for seq_id in missing_queries[:max_show]:

            print(
                f"  missing QUERY: {seq_id}",
                file=sys.stderr
            )

        if len(missing_queries) > max_show:

            print(
                f"  ... and "
                f"{len(missing_queries) - max_show} more",
                file=sys.stderr
            )

    # =============================================================
    # DUPLICATE QUERY WARNING
    # =============================================================

    if duplicate_queries:

        print(
            "",
            file=sys.stderr
        )

        print(
            "WARNING: Some QUERY IDs occur more than once "
            "in the identification table.",
            file=sys.stderr
        )

        print(
            "The FIRST occurrence was used.",
            file=sys.stderr
        )

        max_show = 20

        duplicate_items = list(
            duplicate_queries.items()
        )

        for query, count in duplicate_items[:max_show]:

            print(
                f"  duplicate QUERY: "
                f"{query} ({count} rows)",
                file=sys.stderr
            )

        if len(duplicate_queries) > max_show:

            print(
                f"  ... and "
                f"{len(duplicate_queries) - max_show} more",
                file=sys.stderr
            )

    if (
        not missing_queries
        and
        not duplicate_queries
    ):

        print(
            "",
            file=sys.stderr
        )

        print(
            "No missing or duplicated QUERY IDs detected.",
            file=sys.stderr
        )

    print(
        "",
        file=sys.stderr
    )

    print(
        "Finished.",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()