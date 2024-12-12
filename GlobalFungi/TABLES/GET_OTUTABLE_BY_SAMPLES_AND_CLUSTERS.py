__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import csv
from cStringIO import StringIO  # Python 2

sh_tab = sys.argv[1]            # GF5_RAW_TABLE_TAB_SH.txt.gz
sample_list = sys.argv[2]       # test_samples.txt OR -
shs_list = sys.argv[3]          # test_shs.txt OR -
out_tab = sys.argv[4]           # TAB delimited output OTUTAB

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


#############################################
# LOAD SAMPLES
#############################################
#print("Sample file: " + sample_list)

samples = {}
allSamples = False
i = 0
if sample_list == "-":
    allSamples = True
    print("Samples all")
else:
    for line in openfile(sample_list):
        if not line.rstrip() == "":
            samples[line.rstrip()] = i
            i += 1

print("Samples loaded - " + str(len(samples)))


#############################################
# LOAD SHs
#############################################
#print("SHs file: " + shs_list)

shs = {}
allSHs = False
i = 0
if shs_list == "-":
    allSHs = True
    print("SHs all")
else:
    for line in openfile(shs_list):
        if not line.rstrip() == "":
            shs[line.rstrip()] = i
            i += 1
    print("SHs loaded - " + str(len(shs)))


#############################################
# EXTRACT ABUNDANCES
#############################################
tab_data = {}
h_list = []
s_set = set()
for line in openfile(sh_tab):
        l = line.rstrip()
        vals = l.split('\t')
        sh = vals[0]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        # precheck
        if allSHs or shs.has_key(sh):
            h_list.append(sh)
            abund = {}
            for x in range(len(sn)):
                sample = sn[x]
                if allSamples or samples.has_key(sample):
                    s_set.add(sample)
                    ab = int(sa[x])
                    abund[sample] = ab
            tab_data[sh] = abund

#############################################
# SAMPLES TO LIST
#############################################

s_list = []
for sample in s_set:
    s_list.append(sample)

#############################################
# SAVE OUTPUT
#############################################
# save it
fp = open(out_tab, "w")

# header
line = "SH_name"
for sample in s_list:
    line += "\t" + sample
fp.write(line + "\n")

# body
for sh in h_list:
    line = sh
    abund = tab_data[sh]
    for sample in s_list:
        if abund.has_key(sample):
            line += "\t" + str(abund[sample])
        else:
            line += "\t0"
    fp.write(line + "\n")
 
fp.close()

print("DONE :]")

