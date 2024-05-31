__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import csv
from cStringIO import StringIO  # Python 2

mapping_tab = sys.argv[1]       # NAKI_castleparks_vlk_corrected_derep.map
processed_tab = sys.argv[2]     # NAKI_castleparks_clean_UNITE10_PROCESSED.txt
out_tab = sys.argv[3]           # output
min_sim = float(sys.argv[4])    # similarity threshold - 98.5%
min_cov = float(sys.argv[5])    # coverage threshold - 98%

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
#QUERY   HIT     SIMILARITY      COVERAGE        EVALUE  BITSCORE
#g001239107|size=1       SH0890361.10FU|UDB02107881|k__Fungi;p__Fungi_phy_Incertae_sedis;c__Fungi_cls_Incertae_sedis;o__Fungi_ord_Incertae_sedis;f__Fungi_fam_Incertae_sedis;g__Fungi_gen_Incertae_sedis;s__Fungi_sp|ND**************    95.249  100.0   0.0     656
#g001239106|size=1       SH0890361.10FU|UDB02107881|k__Fungi;p__Fungi_phy_Incertae_sedis;c__Fungi_cls_Incertae_sedis;o__Fungi_ord_Incertae_sedis;f__Fungi_fam_Incertae_sedis;g__Fungi_gen_Incertae_sedis;s__Fungi_sp|ND**************    98.798  100.0   0.0     736


used_variants = {}
taxa_set = set()
for n, line in enumerate(openfile(processed_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        variant = vals[0].split('|')[0]
        sh_name = vals[1].split('|')[0]
        if sh_name != "NO_HIT":
            sim = float(vals[2]) 
            cov = float(vals[3])
            if sim >= min_sim and cov >= min_cov:
                used_variants[variant] = sh_name
                taxa_set.add(sh_name) 

print("Variants passed tresholds - " + str(len(used_variants)))

taxa_list = []
for taxon in taxa_set:
    taxa_list.append(taxon)
taxa_list.append("NO-HIT")

print("Taxons detected - " + str(len(taxa_list)))


#############################################
# compute sample abundance from vars table
#############################################
#g000000001      M03794:84:000000000-CPRCB:1:2111:18070:6901|VJ10
#g000000001      M03794:73:000000000-CJ263:1:2106:24418:11830|VELO04
#g000000001      M03794:84:000000000-CPRCB:1:1111:16036:12903|KK11
#g000000002      M03794:84:000000000-CPRCB:1:1114:7216:15987|VJ12
#g000000002      M03794:84:000000000-CPRCB:1:2105:24183:20757|VJ21
#g000000002      M03794:59:000000000-C57TV:1:1110:6889:24004|JE14
#g000000003      M03794:84:000000000-CPRCB:1:2109:19321:23075|VR15
#g000000004      M03794:48:000000000-BG34K:1:2107:25832:12558|VR15
#g000000005      M03794:48:000000000-BG34K:1:2101:19870:11501|KZ03
#g000000006      M03794:84:000000000-CPRCB:1:1114:18398:21664|KO02


samples = {}
for line in openfile(mapping_tab):
    vals = line.rstrip().split('\t')
    variant = vals[0]
    sample = vals[1].split('|')[1]
    # test if variant (group) pass the thresholds
    taxon = "NO-HIT"
    if used_variants.has_key(variant):
        taxon = used_variants[variant]
    # fill it...
    if samples.has_key(sample):
        taxons = samples[sample]
        if taxons.has_key(taxon):
            taxons[taxon] += 1
        else:
            taxons[taxon] = 1
        samples[sample] = taxons
    else:
        taxons = {}
        taxons[taxon] = 1
        samples[sample] = taxons


print("Variant table processed - samples used " + str(len(samples)))

#############################################
# save the table
#############################################

fp = open(out_tab, "w")
# header
line = "sample"
for taxon in taxa_list:
    line += "\t" + taxon
fp.write(line + "\n")
# body
for sample in samples:
    line = sample
    taxons = samples[sample]
    for taxon in taxa_list:
        if taxons.has_key(taxon):
            line += "\t" + str(taxons[taxon])
        else:
            line += "\t0"
    fp.write(line + "\n")
fp.close()

print("DONE :]")

