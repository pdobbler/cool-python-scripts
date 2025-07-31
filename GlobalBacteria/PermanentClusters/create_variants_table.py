__author__ = 'vetrot'

import sys
import operator
import collections
import gzip
import hashlib

raw_fasta_samples = sys.argv[1] # GB_VOL1_20250526_CLEANED.fa.gz
clustered_vars = sys.argv[2]    # GB_VOL1_20250526_CLEANED_ranked_multi_CLUSTERED.gz - GB_VOL1_20250526_CLEANED_ranked_multi_CLUSTERED_AND_BINNED.gz
clusters_info = sys.argv[3]     # SEEDS_97.0_WORKING_NAMES.fa.info

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def samples_to_strings(samples):
    sample_ids = ";".join(samples.keys())
    counts = ";".join(str(samples[s]) for s in samples.keys())
    return sample_ids, counts

# clustered info
# GB00179323.1      NO_HIT|CL03266|5cbd58e7c6d8d123527eb1ab5719e93b|V_23|S_6|P_1|r_0.00023926054681|SEED
clusters={}
i = 0
for line in openfile(clusters_info, 'r'):
    parts1 = title.split('\t')
    parts2 = parts1.split('|')
    working_name = parts2[0] + '|' + parts2[1]
    clusters[working_name] = parts1[0]
    i += 1

print("Clusters loaded: " + str(i))

# clustered variants
# >p__Verrucomicrobiota|CL050882|8f31d29d77ed4bc76d9a70c9c585f409|V_2|S_1|P_1|r_3.68382469415e-06|97.23320158102767
i = 0
vars_clusters = {}
for line in openfile(clustered_vars, 'r'):
    ch = line[0]
    if ch == '>':
        title = line[1:].strip()
        parts = title.split('|')
        md5_var = parts[2]
        vars_clusters[md5_var] = clusters[parts[0] + '|' + parts[1]]
        i += 1

print("Processed sequence variants: " + str(i))

# raw fasta from samples
# >GB01020442S|An_2019_1acp_Bact|SRR5920425.4|POS=5|POS=253
# TACGTAGGGAGC...
vars_samples = {}
titleRead = False
i = 0
for line in openfile(raw_fasta_samples, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
        i += 1
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            sample_name = title.split('|')[0]
            if vars_samples.has_key(seq):
                samples = vars_samples[seq]
                if samples.has_key(sample_name):
                    samples[sample_name] += 1
                else:
                    samples[sample_name] = 1
                vars_samples[seq] = samples
            else:
                samples = {}
                samples[sample_name] = 1
                vars_samples[seq] = samples

print("Processed raw fasta from samples: " + str(i))



# CREATE TABLE IF NOT EXISTS `variants` (
#   `hash` varchar(32) NOT NULL,
#   `samples` TEXT NOT NULL,
#   `abundances` MEDIUMTEXT NOT NULL,
#   `cluster` varchar(12) NOT NULL,
#   `sequence` TEXT NOT NULL
# );
i = 0
out_file = gzip.open(fasta_file + "_variants.fa.gz", 'wt')
for seq in vars_samples:
    md5_var = hashlib.md5(seq.encode()).hexdigest()
    samples = vars_samples[seq]
    ids_str, counts_str = samples_to_strings(samples)
    out_file.write(md5_var + '\t' + ids_str + '\t' + counts_str + '\t' + vars_clusters[md5_var] + '\t' + seq +'\n')
    i += 1
out_file.close()

print("Variants created: " + str(i))
print("DONE :)")
            





