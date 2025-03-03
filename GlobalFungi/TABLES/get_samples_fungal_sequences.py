__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

annotation = sys.argv[1]    # GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz - info about taxonomy
raw_table = sys.argv[2]    # REL5_REANOT_UNITE10/GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz - info about samples
pair_table = sys.argv[3]    # GF5_RAW_TABLE_SAMPLES.txt.gz

############################################################
# open gzip
############################################################

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def process_variant(samples, sn, sa, sample_pairs):
    for x in range(len(sn)):
        sample = sample_pairs[sn[x]]
        abund = int(sa[x])
        if samples.has_key(sample):
            samples[sample] += abund
        else:
            samples[sample] = abund
    return samples

############################################################
# load pairs id and Sample name
############################################################

# 6043    GF04010997S
pairs = {}
for line in openfile(pair_table, 'r'):
    vals = line.strip().split("\t")
    pairs[vals[0]] = vals[1]


############################################################
# filter taxonomomy to get fungal/nonfungal
############################################################
# QUERY   HIT     SIMILARITY      COVERAGE        EVALUE  BITSCORE
# 3d9468560d2615a3ee0c1a00c1a143b7        SH0953900.10FU  98.718  100.0   6.35e-116       416
# ba029c6222f7dafd1846ff692af82947        SH0980029.10FU  96.622  100.0   1.84e-64        244
# 050932ce68b26bfed574a42f9b66cd85        NO_HIT  -       -       -       -

classification = {}
i = 0
f = 0
for line in openfile(annotation, 'r'):
    if i > 0:
        vals = line.strip().split("\t")
        if vals[1] == 'NO_HIT':
            fungal = False
        else:
            if float(vals[4]) <= 1E-50 or ((float(vals[2])+float(vals[3])) >= 180.0):
                fungal = True
                f += 1
            else:
                fungal = False
        classification[vals[0]] = fungal
    i += 1

print('Annotation processed '+str(i)+" fungal "+str(f))

############################################################
# sample fungal proportion
############################################################
# 4901ef752f778a7001f070ed5362cd1e        56;57;58;59;60;61;62    2;2;11;1;5;1;7  ITS2    0       TAAA...

samples_fungal = {}
samples_nonfun = {}

samples_fungal["ITS1"] = {}
samples_fungal["ITS2"] = {}

samples_nonfun["ITS1"] = {}
samples_nonfun["ITS2"] = {}

for line in openfile(raw_table, 'r'):
    vals = line.strip().split("\t")
    code = vals[0]
    sn = vals[1].split(';') #sample names
    sa = vals[2].split(';') #sample abundances
    marker = vals[3]
    if classification[code]:
        samples_fungal = process_variant(samples_fungal[marker], sn, sa, pairs)
    else:
        samples_nonfun = process_variant(samples_nonfun[marker], sn, sa, pairs)

############################################################
# save output table
############################################################
all_keys = set(samples_fungal.keys()) | set(samples_nonfun.keys())

fp = open(annotation + ".fungalbreakdown", 'w')
fp.write("Sample_ID\tfungal_seqs_ITS1\tnonfungal_seqs_ITS1\ttotal_seqs_ITS1\tfungal_seqs_ITS2\tnonfungal_seqs_ITS2\ttotal_seqs_ITS2\n")
for sample in all_keys:
    # its1
    fun_ITS1 = 0
    if samples_fungal["ITS1"].has_key(sample):
        fun_ITS1 = samples_fungal["ITS1"][sample]
    non_ITS1 = 0
    if samples_nonfun["ITS1"].has_key(sample):
        non_ITS1 = samples_nonfun["ITS1"][sample]

    # its2
    fun_ITS2 = 0
    if samples_fungal["ITS2"].has_key(sample):
        fun_ITS2 = samples_fungal["ITS2"][sample]
    non_ITS1 = 0
    if samples_nonfun["ITS2"].has_key(sample):
        non_ITS2 = samples_nonfun["ITS2"][sample]

    fp.write(sample + "\t" + str(fun_ITS1) + "\t" + str(non_ITS1) + "\t" + str(fun_ITS1 + non_ITS1) + "\t" + str(fun_ITS2) + "\t" + str(non_ITS2) + "\t" + str(fun_ITS2 + non_ITS2) +"\n")
fp.close()

print("DONE :]")
