#!/usr/bin/env python3

import argparse
import gzip
import random
import re
from collections import defaultdict


SIZE_RE = re.compile(r";size=(\d+)")
CLUSTER_RE = re.compile(r"^(CL\d+)\|")


def open_maybe_gzip(filename, mode="rt"):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def parse_fasta(filename):
    header = None
    seq_lines = []

    with open_maybe_gzip(filename, "rt") as f:
        for line in f:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)

                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            yield header, "".join(seq_lines)


def get_cluster_and_size(header):
    cluster_match = CLUSTER_RE.search(header)
    size_match = SIZE_RE.search(header)

    if cluster_match is None:
        raise ValueError(f"Could not find cluster name in header: {header}")

    if size_match is None:
        raise ValueError(f"Could not find size= in header: {header}")

    cluster = cluster_match.group(1)
    size = int(size_match.group(1))

    return cluster, size


def wrap_sequence(seq, width=80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def main():
    parser = argparse.ArgumentParser(
        description="Select one representative FASTA sequence per cluster using largest size="
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file, optionally gzipped"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible tie-breaking"
    )

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    clusters = defaultdict(list)

    for header, seq in parse_fasta(args.input):
        cluster, size = get_cluster_and_size(header)
        clusters[cluster].append({
            "header": header,
            "seq": seq,
            "size": size
        })

    with open(args.output, "w") as out:
        for cluster in sorted(clusters):
            records = clusters[cluster]

            total_size = sum(record["size"] for record in records)
            max_size = max(record["size"] for record in records)

            best_records = [
                record for record in records
                if record["size"] == max_size
            ]

            representative = random.choice(best_records)

            new_header = f">{cluster}|total_size={total_size}|size={representative['size']}"
            out.write(new_header + "\n")

            for line in wrap_sequence(representative["seq"]):
                out.write(line + "\n")


if __name__ == "__main__":
    main()
