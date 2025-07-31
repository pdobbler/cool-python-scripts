__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
from collections import OrderedDict

fasta_file = sys.argv[1]
version = sys.argv[2]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
# >p__Verrucomicrobiota|CL051945|cd51cdb9a4766f7cdda9dee2f5182a99|V_2|S_1|P_1|r_1.49378585086e-05|SEED
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
            rank = float(vals[6].split('_')[1])
            seqs[title] = seq
            ranks[title] = rank

# Keep the sorted dictionary in an ordered way
ordered_ranks = OrderedDict(sorted(ranks.items(), key=lambda x: float(x[1]), reverse=True))

n = 0
digits = 8
fp1 = open(fasta_file+".sorted", 'w')
fp2 = open(fasta_file+".info", 'w')
for title in ordered_ranks:
    n += 1
    name = 'GB';
    for k in range(0,digits-len(str(n))):
        name +='0'
    name += str(n) + '.' + version  
    fp1.write(">" + name + "\n")
    fp1.write(seqs[title] + "\n")
    #    
    fp2.write(name + "\t" + title + "\n")
fp1.close()
fp2.close()

print("Sequences sorted: " + str(n))

