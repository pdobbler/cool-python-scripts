__author__ = 'vetrot'

import sys
import os
import gzip

fasta_file = sys.argv[1]
output_file = sys.argv[2]
# convert_fasta_to_tsv.py

with open(fasta_file, "r") as f, open(output_file, "w") as out:
    seq_id = ""
    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id and seq:
                out.write(f"{seq}\t{seq_id}\n")
            seq_id = line[1:]  # remove ">"
            seq = ""
        else:
            seq += line
    # Write the last record
    if seq_id and seq:
        out.write(f"{seq}\t{seq_id}\n")
