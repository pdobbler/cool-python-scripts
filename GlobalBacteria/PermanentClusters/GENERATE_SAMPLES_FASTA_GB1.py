__author__ = 'vetrot'

import sys
import os
import gzip

variants_table = sys.argv[1] # TABLE_PROCESSED.txt
samples_pairing = sys.argv[2] # GF5_SampleID_to_PermanentID.txt (instead of GF5_RAW_TABLE_SAMPLES.txt.gz)
out_folder = sys.argv[3]

# 660f8b39108cc20d0dbacee2e154f980        16;17   3;3     ITS2    4       4       4       ACACCTC...

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# LOAD PAIRS
#############################################

sample_names = {}
for line in openfile(samples_pairing):
    vals = line.rstrip().split('\t')
    sample_names[vals[0]] = vals[1]
    fp = open(out_folder + vals[1] + "_sample.fasta", 'w')
    fp.close()

print("Sample names loaded... (" + str(len(sample_names)) + ")")

#############################################
# CREATE SEQUENCE FILES
#############################################
n = 0
seqcount = 0
for line in openfile(variants_table):
    vals = line.rstrip().split('\t')
    head = vals[0]
    seq = vals[4]
    samples = vals[1].split(';')
    abunds = vals[2].split(';')
    for i in range(len(samples)):
        sname = sample_names[samples[i]]
        size = abunds[i]
        seqcount += int(size)
        fp = open(out_folder + sname + "_sample.fasta", 'a')
        fp.write(">" + head + "|" + sname + "_size=" + size + "_16S\n")
        fp.write(seq +'\n')
        fp.close()
    n += 1

print "DONE FOR " + str(n) + " VARIANTS - SEQ.COUNT: " + str(seqcount)

#############################################
# CHECK EMPTY FILES / DELETE
#############################################

k = len(sample_names)
for sid in sample_names:
    sname = sample_names[sid]
    sfile = out_folder + sname + "_sample.fasta"
    if os.path.exists(sfile) and os.path.getsize(sfile) == 0:
        # delete empty file...
        os.remove(sfile)
        k -= 1

print "FILLED FILES: " + str(k) + " EMPTY FILES: " + str(len(sample_names)-k) + " (It should be zero!)"




