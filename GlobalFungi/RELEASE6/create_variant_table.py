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

    Returns
    -------
    variants : OrderedDict
        {
            seqID: {
                "sequence": sequence,
                "samples": Counter({sample1: count, sample2: count, ...})
            }
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

        # Sample name = everything before the first "|"
        sample = record_header.split("|", 1)[0]

        # seqID according to the requested method
        seq_id = hashlib.md5(sequence.encode()).hexdigest()

        if seq_id not in variants:
            variants[seq_id] = {
                "sequence": sequence,
                "samples": Counter()
            }
        else:
            # MD5 collision / unexpected inconsistency check
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
                process_record(header, sequence_parts)

                header = line[1:]
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        "FASTA sequence encountered before the first header."
                    )
                sequence_parts.append(line)

        # Last FASTA record
        process_record(header, sequence_parts)

    return variants, fasta_records


def read_identification_table(filename):
    """
    Read gzipped identification table.

    Returns
    -------
    identifications : dict
        QUERY -> {
            "hit": str,
            "similarity": float or None,
            "similarity_text": str,
            "coverage": float or None,
            "coverage_text": str,
            "no_hit": bool
        }

    duplicate_queries : Counter
        QUERY IDs occurring more than once.
    """

    identifications = {}
    query_counts = Counter()

    required_columns = {
        "QUERY",
        "HIT",
        "SIMILARITY",
        "COVERAGE"
    }

    with open_text_gz(filename) as handle:
        header_line = handle.readline()

        if not header_line:
            raise ValueError("Identification table is empty.")

        header = header_line.strip().split()

        missing = required_columns - set(header)

        if missing:
            raise ValueError(
                "Missing required column(s) in identification table: "
                + ", ".join(sorted(missing))
            )

        col = {name: i for i, name in enumerate(header)}

        for line_number, line in enumerate(handle, start=2):
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) < len(header):
                print(
                    f"WARNING: Skipping malformed line {line_number} "
                    f"in identification table.",
                    file=sys.stderr
                )
                continue

            query = fields[col["QUERY"]]
            hit = fields[col["HIT"]]
            similarity_text = fields[col["SIMILARITY"]]
            coverage_text = fields[col["COVERAGE"]]

            query_counts[query] += 1

            # If QUERY occurs more than once, keep only the first row.
            if query in identifications:
                continue

            # ---------------------------------------------------------
            # No HIT case
            # ---------------------------------------------------------
            # Robustly treat a row as "no hit" if HIT, SIMILARITY
            # or COVERAGE contains "-".
            if (
                hit == "NO_HIT"
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

                continue

            # ---------------------------------------------------------
            # Normal HIT
            # ---------------------------------------------------------
            try:
                similarity = float(similarity_text)
                coverage = float(coverage_text)

            except ValueError:
                raise ValueError(
                    f"Invalid SIMILARITY or COVERAGE at line {line_number}: "
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

    return identifications, duplicate_queries


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

    Identification values:
        matching QUERY + thresholds fulfilled -> HIT/similarity/coverage
        matching QUERY + thresholds not fulfilled -> -/-/-
        QUERY absent from identification table -> NA/NA/NA
    """

    if output_filename == "-":
        out = sys.stdout
        close_output = False
    else:
        out = open(output_filename, "w", encoding="utf-8")
        close_output = True

    missing_queries = []

    try:
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

            # Counter preserves insertion order of samples because we
            # populate it while reading the FASTA.
            sample_names = list(data["samples"].keys())

            samples = ";".join(sample_names)

            abundances = ";".join(
                str(data["samples"][sample])
                for sample in sample_names
            )

            if seq_id not in identifications:
                hit = "NA"
                similarity_out = "NA"
                coverage_out = "NA"
                missing_queries.append(seq_id)

            else:
                identification = identifications[seq_id]

                if (
                    identification["similarity"] >= min_similarity
                    and identification["coverage"] >= min_coverage
                ):
                    hit = identification["hit"]
                    similarity_out = identification["similarity_text"]
                    coverage_out = identification["coverage_text"]

                else:
                    hit = "-"
                    similarity_out = "-"
                    coverage_out = "-"

            print(
                "\t".join([
                    seq_id,
                    samples,
                    abundances,
                    hit,
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

    return missing_queries


def main():
    args = parse_arguments()

    print(
        f"Reading FASTA: {args.fasta}",
        file=sys.stderr
    )

    variants, fasta_records = read_fasta(args.fasta)

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

    identifications, duplicate_queries = read_identification_table(
        args.table
    )

    print(
        f"Unique QUERY IDs in identification table: "
        f"{len(identifications)}",
        file=sys.stderr
    )

    missing_queries = write_output(
        variants=variants,
        identifications=identifications,
        min_similarity=args.similarity,
        min_coverage=args.coverage,
        marker=args.marker,
        output_filename=args.output
    )

    # Final report
    print("", file=sys.stderr)
    print("=== FINAL REPORT ===", file=sys.stderr)
    print(
        f"FASTA records:               {fasta_records}",
        file=sys.stderr
    )
    print(
        f"Unique sequence variants:    {len(variants)}",
        file=sys.stderr
    )
    print(
        f"Identification QUERY IDs:    {len(identifications)}",
        file=sys.stderr
    )
    print(
        f"FASTA seqIDs without QUERY:  {len(missing_queries)}",
        file=sys.stderr
    )
    print(
        f"Duplicated QUERY IDs:         {len(duplicate_queries)}",
        file=sys.stderr
    )

    if missing_queries:
        print(
            "\nWARNING: Some FASTA seqIDs were not found as QUERY "
            "in the identification table.",
            file=sys.stderr
        )

        print(
            "Their HIT/SIMILARITY/COVERAGE values were set to NA.",
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
                f"  ... and {len(missing_queries) - max_show} more",
                file=sys.stderr
            )

    if duplicate_queries:
        print(
            "\nWARNING: Some QUERY IDs occur more than once "
            "in the identification table.",
            file=sys.stderr
        )

        print(
            "The FIRST occurrence was used.",
            file=sys.stderr
        )

        max_show = 20

        for query, count in list(duplicate_queries.items())[:max_show]:
            print(
                f"  duplicate QUERY: {query} ({count} rows)",
                file=sys.stderr
            )

        if len(duplicate_queries) > max_show:
            print(
                f"  ... and {len(duplicate_queries) - max_show} more",
                file=sys.stderr
            )

    if not missing_queries and not duplicate_queries:
        print(
            "No missing or duplicated QUERY IDs detected.",
            file=sys.stderr
        )


if __name__ == "__main__":
    main()

