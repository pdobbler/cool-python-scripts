__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import csv
from cStringIO import StringIO  # Python 2

sh_tab = sys.argv[1]
out_tab = sys.argv[2]
sample_list = sys.argv[3]
sample_pairing = sys.argv[4]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
#
#############################################

sample_pairs = {}
for line in openfile(sample_pairing):
     vals = line.rstrip().split('\t')
     sample_pairs[vals[0]] = vals[1]
print("Samples paired - " + str(len(sample_pairs)))


samples_and_shs = {}
for line in openfile(sample_list):
    sample = line.rstrip().split('\t')[0]
    if samples_and_shs.has_key(sample):
        print("duplicite sample: "+sample)
    else:
        samples_and_shs[sample] = {}

print("Samples loaded - " + str(len(samples_and_shs)))


all_shs = set()
for line in openfile(sh_tab):
        l = line.rstrip()
        vals = l.split('\t')
        sh = vals[0]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        #
        for x in range(len(sn)):
            sample = sample_pairs[sn[x]]
            if samples_and_shs.has_key(sample):
                all_shs.add(sh)
                shs = samples_and_shs[sample]
                shs[sh] = int(sa[x])
                samples_and_shs[sample] = shs
                ###

print("Processed found sh: " + str(len(all_shs)))

# change set to list
shs_list = []
for sh in all_shs:
    shs_list.append(sh)


# save it
fp = open(out_tab, "w")
# header
line = "sample"
for sh in shs_list:
    line += "\t" + sh
fp.write(line + "\n")
# body
for sample in samples_and_shs:
    line = sample
    shs = samples_and_shs[sample]
    for sh in shs_list:
        if shs.has_key(sh):
            line += "\t" + str(shs[sh])
        else:
            line += "\t0"
    fp.write(line + "\n")
fp.close()

print("DONE :]")
