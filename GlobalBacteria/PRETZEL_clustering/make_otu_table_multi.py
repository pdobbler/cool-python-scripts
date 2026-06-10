#!/usr/bin/env python3

__author__ = "vetrot"

import argparse
import gzip
import sys
from collections import defaultdict


def openfile(filename, mode="rt"):
    """
    Open plain text or gzipped file.
    Use text mode by default for Python 3.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def fasta_iter(filename):
    """
    Robust FASTA parser.
    Works also if sequences are split across multiple lines.
    Yields: title_without_>, sequence
    """
    title = None
    seq_chunks = []

    with openfile(filename, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if title is not None:
                    yield title, "".join(seq_chunks)

                title = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if title is not None:
            yield title, "".join(seq_chunks)


def load_cluster_variants(cl_vars_files, duplicate_policy="first"):
    """
    Reads one or more clustered variant FASTA files.

    Header example:
    >CL00000001|d829bee4984f82ffc2453212157caf96;samples=25791;relabund_sum=77.2;size=5579362|100.0

    Stores:
    sequence -> CL00000001
    """
    cl_var_seqs = {}
    duplicate_same = 0
    duplicate_conflict = 0

    for cl_vars_file in cl_vars_files:
        print(f"Reading clustered variants: {cl_vars_file}", file=sys.stderr)

        for title, seq in fasta_iter(cl_vars_file):
            cl_name = title.split("|", 1)[0]

            if seq in cl_var_seqs:
                old_cl = cl_var_seqs[seq]

                if old_cl == cl_name:
                    duplicate_same += 1
                    continue

                duplicate_conflict += 1

                if duplicate_policy == "first":
                    continue
                elif duplicate_policy == "last":
                    cl_var_seqs[seq] = cl_name
                elif duplicate_policy == "error":
                    raise ValueError(
                        f"Sequence assigned to multiple clusters: "
                        f"{old_cl} and {cl_name}"
                    )
            else:
                cl_var_seqs[seq] = cl_name

    print(f"Clustered variants loaded: {len(cl_var_seqs)} unique sequences", file=sys.stderr)
    print(f"Duplicate same-cluster variants ignored: {duplicate_same}", file=sys.stderr)
    print(f"Duplicate conflicting variants: {duplicate_conflict}", file=sys.stderr)

    return cl_var_seqs


def count_otus_in_fastas(fasta_files, cl_var_seqs):
    """
    Reads one or more FASTA files with sample names in headers.

    Header example:
    >GB01020442S|An_2019_1acp_Bact|SRR5920425.1|POS=5|POS=253

    sample_name = first field before |
    """
    tab_dict_otus = defaultdict(lambda: defaultdict(int))

    total_records = 0
    matched_records = 0

    for fasta_file in fasta_files:
        print(f"Reading reads FASTA: {fasta_file}", file=sys.stderr)

        for title, seq in fasta_iter(fasta_file):
            total_records += 1

            otu_name = cl_var_seqs.get(seq)
            if otu_name is None:
                continue

            matched_records += 1
            sample_name = title.split("|", 1)[0]
            tab_dict_otus[otu_name][sample_name] += 1

    print(f"Total FASTA records read: {total_records}", file=sys.stderr)
    print(f"Records matched to clusters: {matched_records}", file=sys.stderr)

    return tab_dict_otus


def write_otu_table(tab_dict_otus, output_file):
    with open(output_file, "w") as fp:
        fp.write("OTU\tsample\tabundance\n")

        for otu_name, samples_counts in tab_dict_otus.items():
            samples = []
            abundances = []

            for sample_name, count in samples_counts.items():
                samples.append(sample_name)
                abundances.append(str(count))

            fp.write(
                otu_name
                + "\t"
                + ";".join(samples)
                + "\t"
                + ";".join(abundances)
                + "\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Create compressed OTU table from multiple clustered variant files and multiple FASTA files."
    )

    parser.add_argument(
        "--cl-vars",
        nargs="+",
        required=True,
        help="One or more clustered variant FASTA files, gzipped or plain."
    )

    parser.add_argument(
        "--fasta",
        nargs="+",
        required=True,
        help="One or more original/sample FASTA files, gzipped or plain."
    )

    parser.add_argument(
        "-o",
        "--output",
        default="compressed_otu_tab_otus_to_samples.txt",
        help="Output OTU table."
    )

    parser.add_argument(
        "--duplicate-policy",
        choices=["first", "last", "error"],
        default="first",
        help=(
            "What to do if the same sequence occurs in multiple cluster files "
            "with different cluster names. "
            "'first' keeps the first assignment, 'last' overwrites it, "
            "'error' stops the script."
        )
    )

    args = parser.parse_args()

    cl_var_seqs = load_cluster_variants(
        args.cl_vars,
        duplicate_policy=args.duplicate_policy
    )

    tab_dict_otus = count_otus_in_fastas(
        args.fasta,
        cl_var_seqs
    )

    write_otu_table(tab_dict_otus, args.output)

    print("Done :)", file=sys.stderr)


if __name__ == "__main__":
    main()
    