#!/usr/bin/env python3
"""Create an absolute-count OTU table for selected samples.

Input OTU format (whitespace-separated columns):
    cluster_id    sample1;sample2;...    count1;count2;...

Usage:
    python3 create_absolut_otu_table_sorted_sel_samples.py \
        compressed_otu_table.txt[.gz] list_of_samples.txt output_otu_table.tsv

The output is a tab-separated OTU table with absolute sequence counts.
Only samples listed in ``list_of_samples.txt`` are included. OTUs are sorted
alphabetically by their cluster/OTU identifier. Sample columns retain the order
in which they occur in ``list_of_samples.txt``.
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO


def open_text(filename: str, mode: str = "rt") -> TextIO:
    """Open plain-text or gzip-compressed files based on the .gz suffix."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode, encoding="utf-8")
    return open(filename, mode, encoding="utf-8", newline="")


def read_selected_samples(filename):
    """Read one sample ID per line, preserving order and removing duplicates."""
    samples = []
    seen: set[str] = set()

    with open_text(filename, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            sample = line.strip()

            # Permit empty lines and comment lines in the sample list.
            if not sample or sample.startswith("#"):
                continue

            if sample in seen:
                sys.stderr.write(
                    f"WARNING: duplicate sample {sample!r} on line {line_number} "
                    "of selected-sample file; duplicate ignored\n"
                )
                continue

            seen.add(sample)
            samples.append(sample)

    if not samples:
        raise ValueError(f"No sample IDs found in {filename}")

    return samples


def create_absolute_otu_table(
    otu_file: str,
    selected_samples_file: str,
    output_file: str,
) -> None:
    selected_samples = read_selected_samples(selected_samples_file)
    selected_set = set(selected_samples)
    samples_seen_in_otu: set[str] = set()
    rows = []

    with open_text(otu_file, "rt") as fin:
        for line_number, line in enumerate(fin, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) != 3:
                sys.stderr.write(
                    f"WARNING: line {line_number} skipped; expected 3 columns, "
                    f"got {len(parts)}\n"
                )
                continue

            cluster_id, samples_field, counts_field = parts
            sample_list = samples_field.split(";")
            count_list = counts_field.split(";")

            if len(sample_list) != len(count_list):
                sys.stderr.write(
                    f"WARNING: {cluster_id} on line {line_number} skipped; "
                    "number of samples differs from number of counts\n"
                )
                continue

            absolute_counts: defaultdict[str, int] = defaultdict(int)
            row_is_valid = True

            for sample, count_text in zip(sample_list, count_list):
                if sample not in selected_set:
                    continue

                try:
                    count = int(count_text)
                except ValueError:
                    sys.stderr.write(
                        f"WARNING: {cluster_id} on line {line_number} skipped; "
                        f"invalid integer count {count_text!r} for sample {sample!r}\n"
                    )
                    row_is_valid = False
                    break

                if count < 0:
                    sys.stderr.write(
                        f"WARNING: {cluster_id} on line {line_number} skipped; "
                        f"negative count {count} for sample {sample!r}\n"
                    )
                    row_is_valid = False
                    break

                # Summing also handles an accidental repeated sample in one OTU row.
                absolute_counts[sample] += count
                samples_seen_in_otu.add(sample)

            if not row_is_valid:
                continue

            row = [cluster_id]
            row.extend(absolute_counts[sample] for sample in selected_samples)
            rows.append(row)

    # Keep the same alphabetical OTU ordering as the original script.
    rows.sort(key=lambda row: str(row[0]))

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_text(output_file, "wt") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(["OTU", *selected_samples])
        writer.writerows(rows)

    missing_samples = [
        sample for sample in selected_samples if sample not in samples_seen_in_otu
    ]
    if missing_samples:
        sys.stderr.write(
            "WARNING: selected samples not found in the OTU input: "
            + ", ".join(missing_samples)
            + "\n"
        )

    print(
        f"Done: {output_file} "
        f"({len(rows)} OTUs, {len(selected_samples)} selected samples)"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a tab-separated OTU table with absolute sequence counts "
            "for samples listed in a separate file."
        )
    )
    parser.add_argument(
        "otu_file",
        help="Reduced OTU input: cluster_id, semicolon-separated samples, counts",
    )
    parser.add_argument(
        "sel_samples",
        help="Text file containing one selected sample ID per line",
    )
    parser.add_argument("output_file", help="Output TSV file (optionally .gz)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        create_absolute_otu_table(
            otu_file=args.otu_file,
            selected_samples_file=args.sel_samples,
            output_file=args.output_file,
        )
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
