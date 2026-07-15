#!/usr/bin/env python3

import csv
import sys
import gzip
from collections import defaultdict

otu_file = sys.argv[1]
reads_file = sys.argv[2]
output_file = sys.argv[3]

# Helper: open normal text or .gz text automatically
def open_maybe_gzip(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

# 1) Načtení počtu reads pro každý sample
reads_per_sample = {}

with open(reads_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample = row["sample"]
        reads = int(row["reads"])
        reads_per_sample[sample] = reads

samples_all = sorted(reads_per_sample.keys())

# 2) Převod compressed OTU tabulky na relativní abundance
rows = []

with open_maybe_gzip(otu_file) as fin:
    for line_number, line in enumerate(fin, start=1):
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        if len(parts) != 3:
            sys.stderr.write(
                f"WARNING: line {line_number} skipped, expected 3 columns, got {len(parts)}\n"
            )
            continue

        cluster_id = parts[0]
        sample_list = parts[1].split(";")
        count_list = parts[2].split(";")

        if len(sample_list) != len(count_list):
            sys.stderr.write(
                f"WARNING: {cluster_id} skipped, number of samples != number of counts\n"
            )
            continue

        rel_abund = defaultdict(float)

        for sample, count in zip(sample_list, count_list):
            count = int(count)

            if sample not in reads_per_sample:
                sys.stderr.write(
                    f"WARNING: sample {sample} from {cluster_id} not found in reads_per_sample.tsv\n"
                )
                continue

            total_reads = reads_per_sample[sample]

            if total_reads == 0:
                rel_abund[sample] = 0.0
            else:
                rel_abund[sample] = 100 * count / total_reads

        row = [cluster_id] + [rel_abund[sample] for sample in samples_all]
        rows.append(row)

# Sort rows alphabetically by cluster / OTU name
rows.sort(key=lambda row: row[0])

# Write output
with open(output_file, "w", newline="") as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["OTU"] + samples_all)
    writer.writerows(rows)

print(f"Done: {output_file}")
