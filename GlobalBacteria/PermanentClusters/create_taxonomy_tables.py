__author__ = 'vetrot'

import sys
import os
import gzip
from collections import defaultdict

variant_table = sys.argv[1]         # VARIANTS_TABLE.txt.gz
taxonomy_info = sys.argv[2]         # TAXONOMY_CLUSTERS.txt

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# 11      GB00000011.1    Neobacillus niacini_A_296046    Neobacillus     100.0   100.0   d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Bacillales_B_306087; f__DSM-18226_301387; g__Neobacillus; s__Neobacillus niacini_A_296046 6c5f4f9666d82b1bedf94e99acfe6741
# 60      GB00000060.1    -       -       92.52   100.0   d__Bacteria; p__Chloroflexota; c__Chloroflexia; o__Thermomicrobiales; f__; g__; s__     5b5fc41e4def55055e7e165fecdfc4f1
# 33      GB00000033.1    -       Gp1-AA122       99.209  100.0   d__Bacteria; p__Acidobacteriota; c__Terriglobia; o__Terriglobales; f__SbA1; g__Gp1-AA122; s__   6b315f08e0768a75fac2e33026b2bf5f
cl_sp = {}
cl_gen = {}
cl_names = {}
for line in openfile(taxonomy_info):
    parts = line.rstrip().split("\t")
    cl_names[parts[0]] = parts[1]
    if not parts[2] == '-':
        cl_sp[parts[0]] = parts[2]
    if not parts[3] == '-':
        cl_gen[parts[0]] = parts[3]

print("clusters taxonomy loaded...")

# 2e3427b26a2511095e1ace224422b046        3;1;13;18;2;31;4        1;3;1;2;2;1;1   -       AACG...
# e8043bb6441bce2e4a4b5f8e5c6a6146        37      1       816     TACG...
# cluster_id -> { sample_id: summed_abundance }
cluster_data = defaultdict(lambda: defaultdict(int))
species_data = defaultdict(lambda: defaultdict(int))
genera_data = defaultdict(lambda: defaultdict(int))
cl_var_counts = {}
sp_var_counts = {}
gen_var_counts = {}

for line in openfile(variant_table):
    parts = line.rstrip().split("\t")
    # check structure...
    if len(parts) < 5:
        continue
    # process...
    sample_ids = parts[1].split(';')
    abundances = list(map(int, parts[2].split(';')))
    cluster_id = parts[3]
    # not-binned
    if cluster_id == "-":
        continue
    # clusters
    cl_name = cl_names[cluster_id]
    for sid, ab in zip(sample_ids, abundances):
        cluster_data[cl_name][sid] += ab

    if cl_var_counts.has_key(cl_name):
        cl_var_counts[cl_name] += 1
    else:
        cl_var_counts[cl_name] = 1
    
    # species
    if cl_sp.has_key(cluster_id):
        sp_name = cl_sp[cluster_id]
        for sid, ab in zip(sample_ids, abundances):
            species_data[sp_name][sid] += ab

        if sp_var_counts.has_key(sp_name):
            sp_var_counts[sp_name] += 1
        else:
            sp_var_counts[sp_name] = 1

    # genera
    if cl_gen.has_key(cluster_id):
        gen_name = cl_gen[cluster_id]
        for sid, ab in zip(sample_ids, abundances):
            genera_data[gen_name][sid] += ab

        if gen_var_counts.has_key(gen_name):
            gen_var_counts[gen_name] += 1
        else:
            gen_var_counts[gen_name] = 1

print("variant table processed...")


fp = open("ABUND_TABLE_CLUSTERS.txt", 'w')
for cluster_id, samples in cluster_data.items():
    sample_ids = ';'.join(samples.keys())
    abundances = ';'.join(str(samples[sid]) for sid in samples.keys())
    fp.write(cluster_id + '\t' + sample_ids + '\t' + abundances + '\t' + str(cl_var_counts[cluster_id]) + '\n')
fp.close()

fp = open("ABUND_TABLE_SPECIES.txt", 'w')
for cluster_id, samples in species_data.items():
    sample_ids = ';'.join(samples.keys())
    abundances = ';'.join(str(samples[sid]) for sid in samples.keys())
    fp.write(cluster_id + '\t' + sample_ids + '\t' + abundances + '\t' + str(sp_var_counts[cluster_id]) + '\n')
fp.close()

fp = open("ABUND_TABLE_GENERA.txt", 'w')
for cluster_id, samples in genera_data.items():
    sample_ids = ';'.join(samples.keys())
    abundances = ';'.join(str(samples[sid]) for sid in samples.keys())
    fp.write(cluster_id + '\t' + sample_ids + '\t' + abundances + '\t' + str(gen_var_counts[cluster_id]) + '\n')
fp.close()

print("Done!")