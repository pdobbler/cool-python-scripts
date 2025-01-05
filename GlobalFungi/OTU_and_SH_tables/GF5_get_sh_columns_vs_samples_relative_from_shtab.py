__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

sh_tab = sys.argv[1]                # GF5_RAW_TABLE_TAB_SH.txt.gz
sample_pairs_abund = sys.argv[2]    # GF5_SAMPLES_PAIRS_WITH_TOTAL_ITS.txt
out_tab = sys.argv[3]


#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# CREATE SEQUENCE FILES
#############################################

sample_pairs = {}
sample_abund = {}
for line in openfile(sample_pairs_abund):
     vals = line.rstrip().split('\t')
     sample_pairs[vals[0]] = vals[1]
     sample_abund[vals[0]] = float(vals[2])
print("Samples paired - " + str(len(sample_pairs)))


sample_shs = {}
sh_set = set()
sum = 0
for line in openfile(sh_tab):
        l = line.rstrip()
        vals = l.split('\t')
        sh = vals[0]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        sh_set.add(sh)
        #
        for x in range(len(sn)):
            s_name = sample_pairs[sn[x]]
            s_tot = sample_abund[sn[x]]
            if sample_shs.has_key(s_name):
                shs = sample_shs[s_name]
                shs[sh] = float(sa[x])/s_tot
                sample_shs[s_name] = shs
            else:
                shs = {}
                shs[sh] = float(sa[x])/s_tot
                sample_shs[s_name] = shs

print("sh table processed - total counts " + str(sum) + " num. of samples " + str(len(sh_set)))

# change set to list
sh_list = []
for sh in sh_set:
    sh_list.append(sh)

####################################

fp1 = open(out_tab+".relative.txt", "w")
fp2 = open(out_tab+".presence.txt", "w")
# header
line = "Samples"
for sh in sh_list:
    line += "\t" + sh
fp1.write(line + "\n")
fp2.write(line + "\n")
# body
for s_name in sample_shs:
    line1 = s_name
    line2 = s_name
    shs = sample_shs[s_name]
    for sh in sh_list:
        if shs.has_key(sh):
            line1 += "\t" + str(shs[sh])
            line2 += "\t1"
        else:
            line1 += "\t0"
            line2 += "\t0"
    fp1.write(line1 + "\n")
    fp2.write(line2 + "\n")
# close
fp1.close()
fp2.close()

print("DONE :]")




