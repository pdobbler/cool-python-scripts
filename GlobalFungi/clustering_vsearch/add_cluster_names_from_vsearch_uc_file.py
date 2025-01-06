__author__ = 'vetrot'

import sys
import os

uc_file = sys.argv[1]
fasta_file = sys.argv[2]

i = 0
for line in open(uc_file):
    vals = line.rstrip().split('\t')
    if vals[0] == 'S':
        i += 1
name_max = len(str(i))

print("Clusters count: "+str(i)+" len: "+str(name_max))

sequences = {}
for line in open(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            if sequences.has_key(title):
                print("ERRROR: DUPLICATE TITLE > "+title)
            else:
                sequences[title] = seq

print("Sequences loaded... "+str(len(sequences)))

fpC = open(fasta_file + ".clustered", 'w')
fpS = open(fasta_file + ".seeds", 'w')
for line in open(uc_file):
    vals = line.rstrip().split('\t')
    if vals[0] == 'S':
        cname = vals[1]
        for i in range(name_max - len(cname)):
            cname = "0" + cname
        cname = "PCL" + cname
        #
        title = vals[8]
        # seed
        fpS.write('>' +cname+ '|' + title + '|SEED\n')
        fpS.write(sequences[title] + '\n')
        # renamed seq
        fpC.write('>' +cname+ '|' + title + '|100.0\n')
        fpC.write(sequences[title] + '\n')
    else:
        if vals[0] == 'H':
            cname = vals[1]
            for i in range(name_max - len(cname)):
                cname = "0" + cname
            cname = "PCL" + cname
            #
            title = vals[8]
            sim = vals[3]
            # renamed seq
            fpC.write('>' + cname + '|' + title + '|' + sim + '\n')
            fpC.write(sequences[title] + '\n')
fpC.close()
fpS.close()

print("Saved :)")

