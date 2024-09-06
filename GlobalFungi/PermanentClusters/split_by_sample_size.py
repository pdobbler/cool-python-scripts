__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
from collections import OrderedDict

fasta_file = sys.argv[1]
samples_th = int(sys.argv[2]) # samples threashold

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >fe65866394ee4d161302423df9418134|V_2|S_2|P_2|r_2.0
q = 0
n = 0
fpQ = open(fasta_file+".qualified", 'w')
fpN = open(fasta_file+".nonqualified", 'w')
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
            sample_couts = int(vals[2].split('_')[1])
            if sample_couts >= samples_th:
                fpQ.write(">" + title + "\n")
                fpQ.write(seq + "\n")
                q += 1
            else:
                fpN.write(">" + title + "\n")
                fpN.write(seq + "\n")
                n += 1
fpN.close()
fpQ.close()

print("Sequences saved - qualified variants: "+str(q)+" non-qualified variant: "+str(n))

