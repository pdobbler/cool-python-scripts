__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip


processed_tab = sys.argv[1]     # GF5_RAW_TABLE_PROCESSED_UNITE10.txt
new_output_tab = sys.argv[2]    # output sample to sh

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# read processed tab
#############################################

#CREATE TABLE IF NOT EXISTS `samples_to_sh` (
#  `sample` int NOT NULL,
#  `SHs` MEDIUMTEXT NOT NULL
#);

samples = {}
for line in openfile(processed_tab):
    vals = line.rstrip().split('\t')
    if int(vals[4]) > 0:
        sh_id = vals[4]
        ss = vals[1].split(';')
        for s in ss:
            if samples.has_key(s):
                sh_set = samples[s]
                sh_set.add(sh_id)
                samples[s] = sh_set
            else:
                sh_set = set()
                sh_set.add(sh_id)
                samples[s] = sh_set

print ("samples processed...")

#############################################
# write taxa table
#############################################

fp = open(new_output_tab, "w")

for sample in samples:
    fp.write(sample + "\t" + ';'.join(samples[sample]) + "\n")
fp.close()

print("...SAMPLES TO SH TABLE SAVED")

print("ALL DONE :)")

    