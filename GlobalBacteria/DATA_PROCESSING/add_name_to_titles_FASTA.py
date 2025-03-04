__author__ = 'vetrot'

import sys
import os
import gzip

fasta_file = sys.argv[1]
sample_name = sys.argv[2]

print("name to add to the end of tiitles: " + sample_name)


############################################################
# open gzip
############################################################

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

############################################################


fp = open(fasta_file+".renamed", 'w')
for line in openfile(fasta_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            fp.write(">" + sample_name + '|' + title + "\n")
            fp.write(seq + "\n")
fp.close()


print("Done :]")
