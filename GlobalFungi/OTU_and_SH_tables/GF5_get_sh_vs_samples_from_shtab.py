__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

sh_tab = sys.argv[1]
out_tab = sys.argv[2]
sh_list = sys.argv[3]
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
# CREATE SEQUENCE FILES
#############################################

all_shs = set()
for line in openfile(sh_list):
    vals = line.rstrip().split('\t')
    all_shs.add(vals[0])
print("Shs loaded - " + str(len(all_shs)))


sample_pairs = {}
for line in openfile(sample_pairing):
     vals = line.rstrip().split('\t')
     sample_pairs[vals[0]] = vals[1]
print("Samples paired - " + str(len(sample_pairs)))


shs = {}
samples_set = set()
sum = 0
for line in openfile(sh_tab):
        l = line.rstrip()
        vals = l.split('\t')
        sh = vals[0]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        #
        if sh in all_shs:
            if shs.has_key(sh):
                samples = shs[sh]
                for x in range(len(sn)):
                    sample = sample_pairs[sn[x]]
                    samples_set.add(sample)
                    if samples.has_key(sample):
                        samples[sample] += int(sa[x])
                        sum += int(sa[x])
                    else:
                        samples[sample] = int(sa[x])
                        sum += int(sa[x])
                shs[sh] = samples
            else:
                samples = {}
                for x in range(len(sn)):
                    sample = sample_pairs[sn[x]]
                    samples_set.add(sample)
                    samples[sample] = int(sa[x])
                    sum += int(sa[x])
                shs[sh] = samples

print("sh table processed - total counts " + str(sum) + " num. of samples " + str(len(samples_set)))

# change set to list
samples_list = []
for s in samples_set:
    samples_list.append(s)


sum = 0
# save it
fp = open(out_tab, "w")
# header
line = "SH"
for s in samples_list:
    line += "\t" + s
fp.write(line + "\n")
for sh in shs:
    samples = shs[sh]
    line = sh
    for s in samples_list:
        if samples.has_key(s):
            line += "\t" + str(samples[s])
            sum += samples[s]
        else:
            line += "\t0"
    fp.write(line + "\n")
fp.close()

print("DONE :] Saved counts "+str(sum))



