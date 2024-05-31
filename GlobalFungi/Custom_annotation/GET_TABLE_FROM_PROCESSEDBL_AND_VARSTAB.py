__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import csv
from cStringIO import StringIO  # Python 2

vars_tab = sys.argv[1]          # GF5_RAW_TABLE_PROCESSED.txt.gz
processed_tab = sys.argv[2]     # GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_PROCESSED.txt
out_tab = sys.argv[3]           # output
min_sim = float(sys.argv[4])    # similarity threshold - 98.5%
min_cov = float(sys.argv[5])    # coverage threshold - 98%
sample_pairing = sys.argv[6]      # GF5_RAW_TABLE_SAMPLES.txt.gz

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

used_variants = {}
taxa_set = set()
for n, line in enumerate(openfile(processed_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        sim = float(vals[2]) 
        cov = float(vals[3])
        if sim >= min_sim and cov >= min_cov:
            used_variants[vals[0]] = vals[1]
            taxa_set.add(vals[1]) 

print("Variants passed tresholds - " + str(len(used_variants)))

taxa_list = []
for taxon in taxa_set:
    taxa_list.append(taxon)

print("Taxons detected - " + str(len(taxa_list)))

#############################################
# LOAD PAIRS - sampleID+permanentID 
#############################################
# 75387   GF04000029S

sample_pairs = {}
for line in openfile(sample_pairing):
     vals = line.rstrip().split('\t')
     sample_pairs[vals[0]] = vals[1]

print("Samples paired - " + str(len(sample_pairs)))


#############################################
# compute sample abundance from vars table
#############################################

samples = {}
samples_tot = {}
for line in openfile(vars_tab):
        l = line.rstrip()
        vals = l.split('\t')
        variant = vals[0]
        sn = vals[1].split(';') #sample names
        sa = vals[2].split(';') #sample abundances
        # -----           
        for x in range(len(sn)):
            sample = sample_pairs[sn[x]]
            abund = int(sa[x])
            # sum total abundance
            if samples_tot.has_key(sample):
                samples_tot[sample] += abund
            else:
                samples_tot[sample] = abund
            # check used
            if used_variants.has_key(variant):
                taxon = used_variants[variant]
                if samples.has_key(sample):
                    taxons = samples[sample]
                    if taxons.has_key(taxon):
                        taxons[taxon] += abund
                    else:
                        taxons[taxon] = abund
                    samples[sample] = taxons
                else:
                    taxons = {}
                    taxons[taxon] = abund
                    samples[sample] = taxons


print("Variant table processed - samples total " + str(len(samples_tot)) + " samples used " + str(len(samples)))

#############################################
# save the table
#############################################

fp = open(out_tab, "w")
# header
line = "tot_abund\tsample"
for taxon in taxa_list:
    line += "\t" + taxon
fp.write(line + "\n")
# body
for sample in samples:
    line = str(samples_tot[sample]) + "\t" + sample
    taxons = samples[sample]
    for taxon in taxa_list:
        if taxons.has_key(taxon):
            line += "\t" + str(taxons[taxon])
        else:
            line += "\t0"
    fp.write(line + "\n")
fp.close()

print("DONE :]")

