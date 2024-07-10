__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip


sh_tab = sys.argv[1]            # ELIGIBLE_SAMPLES_ITS1_OUTPUT_TAB.txt
samples_file = sys.argv[2]      # ELIGIBLE_SAMPLES_ITS1.txt
out_tab = sys.argv[3]    # output
abund_thr = float(sys.argv[4])  # 0.00005

# For each SH, samples where its abundance is > 1 (nonsingletons) AND abundances is > 0.00005 are added to the list where SH is present ("1")
# For each SH, samples where its abundance is = 0 are added to the list where SH is absent ("0")
# For each SH, all other samples (i.e., those where the SH was a local singleton or had very low abundance) are excluded from the sample list ("NA"); it is unclear if their low sequence abundance represents a presence or cross-contamination


#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# load samples counts
#############################################

samples_counts = {}
for line in openfile(samples_file):
    vals = line.rstrip().split('\t')
    samples_counts[vals[0]] = int(vals[1])

print("sample counts were loaded... " + str(len(samples_counts)))

#############################################
# compute a new tab!
#############################################
fp = open(out_tab, 'w')

counts_list = []
for n, line in enumerate(openfile(sh_tab)):
    vals = line.rstrip().split('\t')
    if n == 0:
        counts_list.append(0)
        #set samples
        for i in range(1, len(vals)):
            counts_list.append(samples_counts[vals[i]])
    else:
        for i in range(1, len(vals)):
            if int(vals[i]) == 0:
                vals[i] = "0"
            else:
                abund = counts_list[i]/float(vals[i])
                if int(vals[i]) > 1 and abund > abund_thr:
                    vals[i] = "1"
                else:
                    vals[i] = "NA"
    # write...
    fp.write('\t'.join(vals) + "\n")

fp.close()

print("Done :)")
