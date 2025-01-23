__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >CL000012|cf68cb2e631dd3a6c2f730720fec138e|V_13415448|S_11932|P_205|r_67.6273686781|SEED
titleRead = False
fp1 = open(fasta_file+".1", 'w')
fp2 = open(fasta_file+".2", 'w')
for line in openfile(fasta_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            mid = (len(seq) + 1) // 2  # Ensures first half is longer if odd length
            fp1.write(">" + title + "_1\n")
            fp1.write(seq[:mid] + "\n")
            fp2.write(">" + title + "_1\n")
            fp2.write(seq[mid:] + "\n")
fp1.close()
fp2.close()

print("DONE :]")
