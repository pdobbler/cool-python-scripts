__author__ = 'vetrot'

import sys
import os

blast_out = sys.argv[1]

fp = open(blast_out+".shorter.txt", 'w')
for line in open(blast_out):
    parts = line.strip().split("\t")
    parts[1] = parts[1].split("|")[0]
    fp.write("\t".join(parts)+"\n")
fp.close()

print("DONE")


