__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip


original_tab = sys.argv[1]      # GF5_RAW_TABLE_PROCESSED.txt.gz
processed_tab = sys.argv[2]     # GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITEv10_sh_99_s_all_PROCESSED.txt.gz
new_output_tab = sys.argv[3]    # output
min_sim = float(sys.argv[4])    # similarity threshold - 98.5
min_cov = float(sys.argv[5])    # coverage threshold - 90.0

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


#############################################
# load processed blast
#############################################
# 66e0152e23e006fd7e25d2cf3806115e        SH1184039.10FU_MH400205_reps    95.333  100.0   5.86e-62        237

used_variants = {}
for n, line in enumerate(openfile(processed_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        if vals[1] != "NO_HIT":
        	sim = float(vals[2]) 
        	cov = float(vals[3])
	        if sim >= min_sim and cov >= min_cov:
	            used_variants[vals[0]] = vals[1].split('_')[0]

print("Variants passed tresholds - " + str(len(used_variants)))

#############################################
# update processed table
#############################################

sh_map = {}
variants = {}

for n, line in enumerate(openfile(original_tab)):
	vals = line.rstrip().split('\t')
	if used_variants.has_key(vals[0]):
		taxon = used_variants[vals[0]]
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

