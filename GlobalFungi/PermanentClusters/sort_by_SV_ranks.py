__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
from collections import OrderedDict

fasta_file = sys.argv[1]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >fe65866394ee4d161302423df9418134|V_2|S_2|P_2|r_2.0
seqs = {}
ranks = {}
for line in openfile(fasta_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            # process SV
            vals = title.split('|')
            # count samples sequences
            rank = float(vals[4].split('_')[1])
            seqs[title] = seq
            ranks[title] = rank

# Keep the sorted dictionary in an ordered way
ordered_ranks = OrderedDict(sorted(ranks.items(), key=lambda x: float(x[1]), reverse=True))

n = 0
fp = open(fasta_file+".sorted", 'w')
for title in ordered_ranks:
    fp.write(">" + title + "\n")
    fp.write(seqs[title] + "\n")
    n += 1
fp.close()

print("Sequences sorted: " + str(n))

