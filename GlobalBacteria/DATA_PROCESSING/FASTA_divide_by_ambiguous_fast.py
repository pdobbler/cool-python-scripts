__author__ = 'Wietrack 2025'

import sys
import gzip

fasta_file = sys.argv[1]

def openfile(filename, mode='rt'):
    return gzip.open(filename, mode) if filename.endswith('.gz') else open(filename, mode)

first_part = fasta_file[:fasta_file.find('.')]
correct_path = f"{first_part}_correct.fa.gz"
ambiguous_path = f"{first_part}_ambiguous.fa.gz"

valid_bases = {'A', 'C', 'T', 'G'}

i = 0  # correct
n = 0  # ambiguous

with openfile(fasta_file, 'rt') as infile, \
     gzip.open(correct_path, 'wt') as correct_out, \
     gzip.open(ambiguous_path, 'wt') as ambiguous_out:

    while True:
        title = infile.readline()
        seq = infile.readline()
        if not title or not seq:
            break

        title = title.strip()
        seq = seq.strip().upper()

        if set(seq).issubset(valid_bases):
            correct_out.write(f"{title}\n{seq}\n")
            i += 1
        else:
            ambiguous_out.write(f"{title}\n{seq}\n")
            n += 1

total = i + n
print(f"{first_part}{total} sequences processed by ambiguosity - {i} correct sequences vs {n} ambiguous sequences")
