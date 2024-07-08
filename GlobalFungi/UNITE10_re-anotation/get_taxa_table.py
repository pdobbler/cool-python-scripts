__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip


processed_tab = sys.argv[1]     # GF5_RAW_TABLE_PROCESSED_UNITE10.txt
taxonomy_tab = sys.argv[2]      # UPDATED_TAX_TABLE.txt
new_output_tab = sys.argv[3]    # output
taxa_level = int(sys.argv[4])   # taxa level - 7-species, 6-genus

print("Taxa laevel: "+str(taxa_level))

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# load taxonomy tab
#############################################

level_pair = {}
for line in openfile(taxonomy_tab):
    vals = line.rstrip().split('\t')
    level_pair[vals[8]] = vals[taxa_level]

#############################################
# read processed tab
#############################################

sh_map = {}
variants = {}
for line in openfile(processed_tab):
    vals = line.rstrip().split('\t')
    if int(vals[4]) > 0:
        taxon = level_pair[vals[4]]
        ss = vals[1].split(';')
        aa = vals[2].split(';')    
        if sh_map.has_key(taxon):
            samples = sh_map[taxon]
            for i in range(len(ss)):
                if samples.has_key(ss[i]):
                    samples[ss[i]] += int(aa[i])
                else:
                    samples[ss[i]] = int(aa[i])
            sh_map[taxon] = samples
            variants[taxon] += 1
        else:
            samples = {}
            for i in range(len(ss)):
                samples[ss[i]] = int(aa[i])
            sh_map[taxon] = samples
            variants[taxon] = 1

print ("taxa processed...")

#############################################
# write taxa table
#############################################

fp = open(new_output_tab, "w")

for taxon in sh_map:
    samples = sh_map[taxon]
    ss = []
    aa = []
    for s_id in samples:
        ss.append(s_id)
        aa.append(str(samples[s_id])) 
    line = taxon + "\t" + ';'.join(ss) + "\t" + ';'.join(aa) + "\t" + str(variants[taxon])
    fp.write(line + "\n")
fp.close()

print("...TAXA TABLE SAVED")

print("ALL DONE :)")

    