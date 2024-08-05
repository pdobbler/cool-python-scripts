__author__ = 'vetrot'

import sys
import os

fasta_file = sys.argv[1]

names = set()
variants = {}
sequences = {}
i = 0
for line in open(fasta_file):
    ch = line[0]
    if ch == '>':
        i += 1
        titleRead = True
        title = line[1:].strip()
        name = title.split("|")[0]
        if variants.has_key(name):
            variants[name] += 1
            sequences[name] += int(title.split("|")[2].split("_")[1])
        else:
            names.add(name)
            variants[name] = 1
            sequences[name] = int(title.split("|")[2].split("_")[1])

print("Sequences loaded "+str(i)+" clusters found "+str(len(names)))

fp = open(fasta_file+".counts.txt", 'w')
fp.write('CLUSTER_NAME\tVAR_SIZE\tSEQ_SIZE\n')
for name in names:
    fp.write(name+'\t'+str(variants[name])+'\t'+str(sequences[name])+'\n')
fp.close()

print("DONE")