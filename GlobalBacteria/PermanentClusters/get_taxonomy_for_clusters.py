__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

processed_taxonomy = sys.argv[1]    # GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt.gz
taxonomy_tab = sys.argv[2]          # taxonomy.tsv
seeds_fasta = sys.argv[3]           # SEEDS_97.0_WORKING_NAMES.fa
composed_out = sys.argv[4]          # output


#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


#############################################
# LOAD TAXONOMY
#############################################
taxa = {}
for line in openfile(taxonomy_tab):
    parts = line.rstrip().split('\t')
    taxa[parts[0]] = parts[1]

print("Taxonomy loaded...")

#############################################
# LOAD SEEDs
#############################################
# >p__Pseudomonadota|CL000005|6e105f3f79341d4c6ac024b928786898|V_1211301|S_8765|P_81|r_12.7105402863|SEED
# ACTG...
seeds = {}
for line in openfile(seeds_fasta):
    if line[0] == '>':
        parts = line[1:].rstrip().split('|')
        seeds[parts[2]] = parts[1] + '_' + parts[0] + '\t' + parts[3].split('_')[1] + '\t' + parts[4].split('_')[1] + '\t' + parts[5].split('_')[1] + '\t' + parts[6].split('_')[1]

print("SEEDs loaded...")

#############################################
# CHECK TAXONOMY
#############################################

# Uniq13;size=1651634     f61db4f18408220f35b368effea7821f        G003167215      95.652  100.0   2.61e-112       407     UNKNOWN 16S     UNKNOWN
fp = open(composed_out, "w")
fp.write("variant\tcluster\tvar_size\tvar_samples\tvar_papers\tvar_rank\tsimilarity\tcoverage\tfull_taxonomy\n")
for line in openfile(processed_taxonomy):
    parts = line.rstrip().split('\t')
    if seeds.has_key(parts[1]):
        if parts[2] == "NO_HIT":
            fp.write(parts[1] + "\t"+ seeds[parts[1]] + "\t-\t-\tNO_HIT\n")
        else:
            fp.write(parts[1] + "\t"+ seeds[parts[1]] + "\t" + parts[3] + "\t" + parts[4] + "\t" + taxa[parts[2]] + "\n")
fp.close()

print("Done :)")



