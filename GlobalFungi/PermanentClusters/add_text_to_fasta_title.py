__author__ = 'vetrot'

import sys
import os
import gzip

fasta_file = sys.argv[1]
text = sys.argv[2]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

print("name to add to start of titles: " + sample_name)

i = 0
filled = False
fp = open(fasta_file + ".renamed", 'w')
for n, line in enumerate(openfile(fasta_file)):
    if n % 2 == 0:
        r1_0 = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            r1_1 = line.rstrip()
            filled = True
    if filled:
        fp.write(">" + text + r1_0 + '\n')
        fp.write(r1_1 + '\n')
        filled = False
fp.close()

print("Done :]")
