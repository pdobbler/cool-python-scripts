__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

orig_fasta = sys.argv[1]    # GF5_ITS2_PERM_CLUSTERS_SEEDS_ALL_SORTED_RENAMED.fa
chim_fasta = sys.argv[2]    # chimeras.fasta
out_file = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >PC00000029|092b793e5c586b9fe88f8e7585cdf765|V_438969|S_1795|P_80|r_39.898434862
chimeras = {}
for line in openfile(chim_fasta, 'r'):
    ch = line[0]
    if ch == '>':
        title = line[1:].strip()
        chimeras[title] = "chimeric"


# >PC00000012|cf68cb2e631dd3a6c2f730720fec138e|V_13415448|S_11932|P_205|r_67.6273686781
fp = open(out_file, 'w')
for line in openfile(orig_fasta, 'r'):
    ch = line[0]
    if ch == '>':
        title = line[1:].strip()
        state = "OK"
        if chimeras.haskey(title):
            state = "chimeric"
        fp.write(title + "\t" + state + "\n")
fp.close()


print("Chimeric status created...")


