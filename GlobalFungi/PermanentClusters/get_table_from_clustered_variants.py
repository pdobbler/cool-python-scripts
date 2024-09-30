__author__ = 'avetrot'

import sys
import gzip
import hashlib

CLUSTERED_VARS = sys.argv[1] #REL4_ITS2_FINAL_qualified_clustered_and_binned.fa
RAW_GF_FASTA = sys.argv[2]   #REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >PCL000001|ad765e6d39aa1c7a3811388dd6c850c6|V_33910054|S_12224|P_146|r_0.728005994798|100.0

variants = {}
clusters = {}
for line in openfile(CLUSTERED_VARS, 'r'):
    ch = line[0]
    if ch == '>':
        parts = line[1:].strip().split('|')
        variants[parts[1]] = parts[0]
        if not clusters.has_key(parts[0]):
            clusters[parts[0]] = {}

print("Sequence variants loaded... "+str(len(variants)))

# >AMF007119v1|Anthony_2020_mBY|S1_
# >GF4S07047b|Sun_2021_PK|ERR4887914.1283241

titleRead = False
samples = set()
for line in openfile(RAW_GF_FASTA, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            md5 = hashlib.md5(line.strip().encode()).hexdigest()
            if variants.has_key(md5):
                sample_name = title.split('|')[0]
                samples.add(sample_name)
                samples = clusters[variants[md5]]
                if samples.has_key(sample_name):
                    samples[sample_name] += 1
                else:
                    samples[sample_name] = 1
                clusters[variants[md5]] = samples

# Convert the set to a list
samples_list = list(samples)

# save output
fp = open(CLUSTERED_VARS+".TABLE.txt", 'w')
line = "cluster_name"
for sample_name in samples_list:
    line += "\t" + sample_name
fp.write(line + '\n')

for cl_name in clusters:
    line = "cluster_name"
    for sample_name in samples_list:
        samples = clusters[sample_name]
        if samples.has_key(sample_name):
            line += "\t" + str(samples[sample_name])
        else:
            line += "\t0"
    fp.write(line + '\n')
fp.close()

print("Table saved...")




