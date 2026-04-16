__author__ = 'Wietrack 2019, modified 2026'

import sys
import hashlib
import gzip
from collections import defaultdict

fasta_file = sys.argv[1]

################################################
def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
################################################

# total number of reads per sample
sample_total_reads = defaultdict(int)

# counts of each sequence in each sample:
# seq_sample_counts[seq][sample] = count
seq_sample_counts = defaultdict(lambda: defaultdict(int))

i = 0
current_sample = None
title_read = False

with openfile(fasta_file, 'rt') as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue

        if line.startswith('>'):
            title_read = True

            # header example: >GB01001080S|seq9999
            header = line[1:]
            sample_name = header.split('|', 1)[0]
            current_sample = sample_name
        else:
            if title_read:
                i += 1
                title_read = False

                seq = line
                sample_total_reads[current_sample] += 1
                seq_sample_counts[seq][current_sample] += 1

print("Processed sequences " + str(i) + " to variants " + str(len(seq_sample_counts)))

# Build stats for each sequence variant
# item = (seq, sample_count, relabund_sum, total_size)
variant_stats = []

for seq, per_sample_counts in seq_sample_counts.items():
    total_size = sum(per_sample_counts.values())
    sample_count = len(per_sample_counts)

    relabund_sum = 0.0
    for sample, count in per_sample_counts.items():
        relabund_sum += float(count) / float(sample_total_reads[sample])

    variant_stats.append((seq, sample_count, relabund_sum, total_size))

# Sort:
# 1) by number of samples descending
# 2) by sum of relative abundances descending
# 3) optionally by total size descending for deterministic output
variant_stats.sort(key=lambda x: (x[1], x[2], x[3]), reverse=True)

print("Variants sorted...")

md5_titles = {}
out_file = gzip.open(fasta_file + "_scored_variants.fa.gz", 'wt')

for seq, sample_count, relabund_sum, total_size in variant_stats:
    md5_title = hashlib.md5(seq.encode()).hexdigest()

    out_file.write(
        ">" + md5_title +
        ";samples=" + str(sample_count) +
        ";relabund_sum=" + "{:.10f}".format(relabund_sum) +
        ";size=" + str(total_size) + "\n"
    )
    out_file.write(seq + "\n")

    if md5_title in md5_titles:
        print("Error: nonunique md5 " + md5_title + " seq: " + seq + " size " + str(total_size) + " !!!")
        md5_titles[md5_title] += 1
    else:
        md5_titles[md5_title] = 1

out_file.close()

print("Done :)")
