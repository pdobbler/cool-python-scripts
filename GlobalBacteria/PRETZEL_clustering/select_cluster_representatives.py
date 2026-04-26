#!/usr/bin/env python3

import argparse
import gzip
import random
import re
from collections import defaultdict


SIZE_RE = re.compile(r";size=(\d+)")
CLUSTER_RE = re.compile(r"^(CL\d+)\|")
MD5_RE = re.compile(r"^CL\d+\|([^;|]+)")


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


def get_cluster_md5_and_size(header):
    cluster_match = CLUSTER_RE.search(header)
    md5_match = MD5_RE.search(header)
    size_match = SIZE_RE.search(header)

    if cluster_match is None:
        raise ValueError(f"Could not find cluster name in header: {header}")

    if md5_match is None:
        raise ValueError(f"Could not find MD5/hash in header: {header}")

    if size_match is None:
        raise ValueError(f"Could not find size= in header: {header}")

    cluster = cluster_match.group(1)
    md5 = md5_match.group(1)
    size = int(size_match.group(1))

    return cluster, md5, size


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
        help="Output FASTA file. Use .gz extension to write gzipped output."
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
        cluster, md5, size = get_cluster_md5_and_size(header)

        clusters[cluster].append({
            "seq": seq,
            "md5": md5,
            "size": size
        })

    with open_maybe_gzip(args.output, "wt") as out:
        for cluster in sorted(clusters):
            records = clusters[cluster]

            total_size = sum(record["size"] for record in records)
            max_size = max(record["size"] for record in records)

            best_records = [
                record for record in records
                if record["size"] == max_size
            ]

            representative = random.choice(best_records)

            new_header = (
                f">{cluster}"
                f"|{representative['md5']}"
                f"|total_size={total_size}"
                f"|representative_size={representative['size']}"
            )

            out.write(new_header + "\n")
            out.write(representative["seq"] + "\n")


if __name__ == "__main__":
    main()
