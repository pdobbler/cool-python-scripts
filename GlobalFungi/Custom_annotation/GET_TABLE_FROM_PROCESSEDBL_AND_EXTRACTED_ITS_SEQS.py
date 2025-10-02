__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

extracted_fasta = sys.argv[1]   # GSSP_CLEANED_EXTRACTED_ITS2.fasta.gz
processed_tab = sys.argv[2]     # GSSP_CLEANED_EXTRACTED_ITS2_uniq_UNITE10ECO_PROCESSED.txt
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
# load fasta
#############################################


#############################################
# load processed blast
#############################################

# SEQ_NAME        SEQ_HASH        HIT     SIMILARITY      COVERAGE        EVALUE  BITSCORE        SH      MARKER  KINGDOM
# Uniq1;size=10067879     82780436adb026d0faf0ab5680b15b13        SH0962330.10FU|EF679363|k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Cladosporiales;f__Cladosporiaceae;g__Cladosporium;s__Cladosporium_herbarum|litter_saprotroph*plant_pathogen*foliar_endophyte*leaf/fruit/seed_pathogen*wood,leaf/fruit/seed,soil,dung,animal_material**partly_aquatic*animal_parasite*filamentous_mycelium*perithecium_(hymenium_hidden,_narrow_opening)*closed**** 100.000 100.0   2.39e-73        274     SH0962330.10FU  ITS2    k__Fungi
# Uniq2;size=7967525      e399dc7c0823811b06f9eea10fa4ea6c        SH0962330.10FU|EF679363|k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Cladosporiales;f__Cladosporiaceae;g__Cladosporium;s__Cladosporium_herbarum|litter_saprotroph*plant_pathogen*foliar_endophyte*leaf/fruit/seed_pathogen*wood,leaf/fruit/seed,soil,dung,animal_material**partly_aquatic*animal_parasite*filamentous_mycelium*perithecium_(hymenium_hidden,_narrow_opening)*closed**** 99.324  100.0   1.11e-71        268     SH0962330.10FU  ITS2    k__Fungi
# Uniq519921;size=8       115196452979d6e3d74a9a4aefa75a90        NO_HIT  -       -       -       -       NO_HIT  ITS2    NO_HIT

variant_tax = {}
for n, line in enumerate(openfile(processed_tab)):
    if n>0: 
        vals = line.rstrip().split('\t')
        if vals[2] != 'NO_HIT':
            sim = float(vals[3]) 
            cov = float(vals[4])
            if sim >= min_sim and cov >= min_cov:
                variant_tax[vals[1]] = vals[2].split('|')[0]

print("blast variants loaded...")

# >ERR12552472.136614|POS=1|POS=316
# CACCACTCAAGCCTCGCTTGGTATTGGGCAACGCGGTCCGCCGCGTGCCTCAAATCGTCCGGCTGGGTCTTCTGTCCCCTAAGCGTTGTGGAAACTATTCGCTAAAGGGTGTTCGGGAGGCTACGCCGTAAAACAACCCCATTTCTAA
# >ERR12552544.158152|POS=0|POS=322
# CACCACTCAAGCTATGCTTGGTATTGGGCGTCGTCCTTAGTTGGGCGCGCCTTAAAGACCTCGGCGAGGCCACTCCGGCTTTAGGCGTAGTAGAATTTATTCGAACGTCTGTCAAAGGAGAGGAACTCTGCCGACTGAAACCTTTATTTTTCTA

sh_samples = {}
sample_names = set()
titleRead = False
for line in openfile(extracted_fasta):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            gname = hashlib.md5(seq.encode()).hexdigest()
            if variant_tax.has_key(gname):
                sample = title.split('.')[0]
                sample_names.add(sample)
                sh_name = variant_tax[gname]
                if sh_samples.has_key(sh_name):
                    samples = sh_samples[sh_name]
                    if samples.has_key(sample):
                        samples[sample] += 1
                    else:
                        samples[sample] = 1
                    sh_samples[sh_name] = samples
                else:
                    samples = {}
                    samples[sample] = 1
                    sh_samples[sh_name] = samples

print("FASTA processed - num of SH: " + str(len(sh_samples)) + " num of samples: " + str(len(sample_names)))

sample_list = []
for sample in sample_names:
    sample_list.append(sample)

print("Samples found: " + str(len(sample_list)))


#############################################
# save otutable
#############################################

fp = open(out_tab, "w")
line = "SH"
for sample in sample_list:
    line += "\t" + sample
fp.write(line + '\n')
for sh_name in sh_samples:
    samples = sh_samples[sh_name]
    line = sh_name
    for sample in sample_list:
        if samples.has_key(sample):
            line += "\t" + str(samples[sample])
        else:
            line += "\t0"
    fp.write(line + '\n')
fp.close()

print("DONE :]")

