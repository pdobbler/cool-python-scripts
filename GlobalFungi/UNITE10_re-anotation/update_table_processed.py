__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip


original_tab = sys.argv[1]      # GF5_RAW_TABLE_PROCESSED.txt.gz
processed_tab = sys.argv[2]     # GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz
new_output_tab = sys.argv[3]    # output
min_sim = float(sys.argv[4])    # similarity threshold - 98.5
min_cov = float(sys.argv[5])    # coverage threshold - 90.0
complete_tax_tab = sys.argv[6]	# unite10_complete_taxonomy.txt

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# load complete taxonomy
#############################################

complete_taxonomy = {}
for n, line in enumerate(openfile(complete_tax_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        complete_taxonomy[vals[0]] = vals[1] + '|' + vals[2] + '|' + vals[3] + '|' + vals[4] + '|' + vals[5] + '|' + vals[6] + '|' + vals[7]


print("Complete taxonomy loaded - " + str(len(complete_taxonomy)))

#############################################
# load processed blast
#############################################

used_variants = {}
for n, line in enumerate(openfile(processed_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        if vals[1] != "NO_HIT":
        	sim = float(vals[2]) 
        	cov = float(vals[3])
	        if sim >= min_sim and cov >= min_cov:
	            used_variants[vals[0]] = vals[1]

print("Variants passed tresholds - " + str(len(used_variants)))

#############################################
# update processed table
#############################################

# write output table
fp = open(new_output_tab, "w")

sh_map = {}
for n, line in enumerate(openfile(original_tab)):
	vals = line.rstrip().split('\t')
	line = vals[0] + "\t" + vals[1] + "\t" + vals[2] + "\t" + vals[3]
	if used_variants.has_key(vals[0]):
		sh = used_variants[vals[0]]
		if sh_map.has_key(sh):
			sh_id = sh_map[sh]
		else:
			sh_id = len(sh_map) + 1
			sh_map[sh] = sh_id
		line += "\t" + str(sh_id)
	else:
		line += "\t0"
	line += "\t" + vals[7]
	fp.write(line + "\n")

fp.close()
print("...ROCESS TABLE UPDATED")

#############################################
# write taxonomy table
#############################################

fp = open("UPDATED_TAX_TABLE.txt", "w")

# Sort the dictionary by value
sorted_sh_map = dict(sorted(sh_map.items(), key=lambda item: item[1]))

# Iterate over the sorted dictionary
for key, value in sorted_sh_map.items():
    #print(f"{key}: {value}")
    taxonomy = complete_taxonomy[key].split('|')
    line = key + '\t' + '\t'.join(taxonomy) + '\t' + str(value)
    fp.write(line + "\n")
fp.close()
print("...TAXONOMY EXPORTED")

print("ALL DONE :)")



