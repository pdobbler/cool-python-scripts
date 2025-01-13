__author__ = 'avetrot'

import sys
import hashlib

FASTA = sys.argv[1]
PROCESSED_BLAST = sys.argv[2]
OUT_TABLE = sys.argv[3]

binned={}
for n, line in enumerate(open(PROCESSED_BLAST)):
    if n>0:
        vars = line.rstrip().split('\t')
        binned[vars[1]] = vars[2]

print("hashcodes and VTX were set...")

table = {}
filled = False
vtx_names = {}
i = 0
k = 0
for n, line in enumerate(open(FASTA)):
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line.rstrip()
            filled = True
    if filled:
        filled = False
        sample = title.split('|')[0]
        hash_name = hashlib.md5(seq.encode()).hexdigest()
        if  binned.has_key(hash_name):
            i += 1
            vtx = binned[hash_name]
            vtx_names[vtx] = 0
            if table.has_key(sample):
                vtxs = table[sample]
                if vtxs.has_key(vtx):
                    vtxs[vtx] = vtxs[vtx] + 1
                else:
                    vtxs[vtx] = 1
                table[sample] = vtxs
            else:
                vtxs = {vtx: 1}
                table[sample] = vtxs
        else:
            k = k + 1

print("FASTA was processed - total seqs: "+str(k+i)+" - classified "+str(i)+" - not classified "+str(k)+" - #vtx: "+str(len(vtx_names)))

fp = open(OUT_TABLE, 'w')
# header
line = "samples"
for vtx in vtx_names:
    line += "\t" + vtx
fp.write(line + '\n')
for sample in table:
    line = sample
    vtxs = table[sample]
    for vtx in vtx_names:
        if vtxs.has_key(vtx):
            line += "\t" + str(vtxs[vtx])
        else:
            line += "\t0"
    fp.write(line + '\n')
fp.close()

print("DONE :)")



