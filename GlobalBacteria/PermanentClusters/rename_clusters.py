__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import re

cluster_names = sys.argv[1]   # /mnt/DATA1/GLOBAL_BACTERIA/FINAL/SEEDS_97.0_WORKING_NAMES.fa.info
input_fasta = sys.argv[2]     # GB_VOL1_20250526_CLEANED_ranked_multi_CLUSTERED.gz
out_fasta = sys.argv[3] 

#############################################
# HELPERS
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# LOAD CLUSTER NAMES
#############################################

# GB00000001.1    p__Pseudomonadota|CL000001|6472eb8b1e09f892aca2f23182962903|V_4505598|S_15048|P_97|r_65.6430231571|SEED
# GB00000002.1    p__Actinomycetota|CL000001|7ff346973a282aa55de296afdb5d74af|V_3348910|S_11165|P_91|r_46.9649270943|SEED

clusters = {}
for line in openfile(cluster_names):
    vals = line.rstrip().split('\t')
    parts = vals[1].split('|')
    clusters[parts[0] + "|" + parts[1]] = vals[0]

print("Cluster names loaded...")

# >NO_HIT|CL00002|76d6776fc11f3066dbab47575a4c7c6f|V_12826|S_154|P_4|r_0.0336183871779|97.11934156378601
# TACGTGAGAGGCAAGCGTTATTCGTCATTAATGGGTCTAAAGGGTACGTAGGCGGTATAGTAAATCTTTTCTTAAATAACACTTAAAAGTGGGTTTGATAATGCTAAACTAGAGTTAAAAGGAGTAACAAATACGACAAGTGGAGTGTAACAATACTTAGATACTTGAAGGGTTGCGAAGTGGCGAAGGCATGTTACTATTGAAAACTGACGCTGAGGTACGAAGGCATGGGTATCGATCGGG
# >NO_HIT|CL00010|df6eca6d40595fa261bef410b75af8fb|V_2436|S_83|P_8|r_0.0332094675183|97.93388429752066
# TACGTGAGAGGCAAGCGTTATTCGTCATTAATGGGTCTAAAGGGTACGTAGGCGGTATAAAAAATCTTATTTGTAGACGTGAGAGGTGAGTAATGATAAATTTATACTAGAGTCAGAAAGGAGTAACAAATACAACGAGAGGAGTGTAACAATACGTAGATACTCGAGAGGTTGCGAAGTGGCGAAGGCGTGTTACTATTGATGACTGACGCTGAGGTACGAAGGCATGGGTATCGAACGGG

of = gzip.open(out_fasta + ".gz", 'wt')
# load fasta seqs
titleRead = False
for line in openfile(input_fasta):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            parts = title.rstrip().split('|', 2)  # split into first two + rest
            prefix = parts[0] + "|" + parts[1]
            of.write('>'+ clusters[prefix] + "|" + parts[2] + "\n")
            of.write(seq + '\n')

print("FASTA loaded...")

of.close()

print("DONE :)")
