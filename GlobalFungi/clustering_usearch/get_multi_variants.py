__author__ = 'avetrot'

import sys
import gzip

FASTA_in = sys.argv[1]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# >Uniq1;size=8939333
# >UniqXXX;size=1

fp = open(FASTA_in+".multi", 'w')
fs = open(FASTA_in+".single", 'w')
variants={}
filled = False
for n, line in enumerate(openfile(FASTA_in)):
    if n % 2 == 0:
        r1_0 = line.rstrip()
    else:
        if n % 2 == 1:
            r1_1 = line.rstrip()
            filled = True
    if filled:
        if r1_0.split(';')[1] == "size=1":
            fs.write(r1_0 + '\n')
            fs.write(r1_1 + '\n')
        else:
            fp.write(r1_0 + '\n')
            fp.write(r1_1 + '\n')
        filled = False
fp.close()
fs.close()
print("Done :)")

