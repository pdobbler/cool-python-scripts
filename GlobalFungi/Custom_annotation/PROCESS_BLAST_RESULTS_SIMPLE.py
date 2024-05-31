__author__ = 'Wietrack 2023'

import sys
import os
import gzip

fasta_file = sys.argv[1]
blast_out6 = sys.argv[2]
out_table = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# load blast results...
blast_info = {}
for line in openfile(blast_out6):
    values = line.strip().split("\t")
    # storing - pident | length |  evalue | bitscore
    blast_info[values[0]] = values[1] + "\t" + values[2] + "\t" + str(int(values[7])-int(values[6])+1) + "\t" + values[10] + "\t" + values[11]

# loop over fasta...
ot = open(out_table, 'w')
loaded = False
titleRead = False
i = 0
n = 0
f = 0
#load fasta sequences
ot.write("QUERY\tHIT\tSIMILARITY\tCOVERAGE\tEVALUE\tBITSCORE\n")
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            slen = len(seq)
            if slen>0:
                if blast_info.has_key(title):
                    values = blast_info[title].split("\t")
                    ot.write(title+"\t"+values[0]+"\t"+values[1]+"\t"+str((float( values[2])/slen)*100)+"\t"+ values[3]+"\t"+ values[4]+"\n")
                    f = f + 1
                else:
                    ot.write(title+"\tNO_HIT\t-\t-\t-\t-\n")
                i=i+1
            else:
                n = n + 1
ot.close()


print(str(i)+" sequences loaded correctly - "+str(n)+" sequnces are empty - Omitted!")
print("FOUND BLAST RESULTS FOR "+str(f)+" SEQUENCES OUT OF "+str(i))
