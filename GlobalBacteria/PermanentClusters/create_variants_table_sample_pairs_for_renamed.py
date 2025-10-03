__author__ = 'vetrot'

import sys
import operator
import collections
import gzip
import hashlib

raw_fasta_samples = sys.argv[1] # GB_VOL1_20250526_CLEANED.fa.gz
clustered_vars = sys.argv[2]    # GB_VOL1_20250526_CLEANED_ranked_ALL_CLUSTERED_RENAMED.gz

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def samples_to_strings(samples):
    sample_ids = ";".join(str(sample_names[k]) for k in samples.keys())
    counts = ";".join(str(v) for v in samples.values())
    return sample_ids, counts

# Input: dictionary of sample occurrences
# Example: sample_occurence = {"S1": 10, "S2": 5, "S3": 15}
def create_sample_names(sample_occurence):
    # Sort sample names by decreasing occurrence count
    sorted_samples = sorted(sample_occurence.items(), key=lambda x: x[1], reverse=True)
    # Create dictionary mapping sample_name -> index (starting from 1)
    sample_names = {}
    for idx, (name, count) in enumerate(sorted_samples, start=1):
        sample_names[name] = idx    
    return sample_names

def hit_numeric_key(hit):
    """
    Extract the main numeric ID from a hit name like 'GB00002469.1' -> 2469.
    If multiple numeric groups are present, ignore the final version suffix.
    """
    import re
    # Match GB + digits + .version
    m = re.match(r'^[A-Z]+0*(\d+)\.\d+$', hit)
    if m:
        return int(m.group(1))
    # Fallback: take first integer group
    m = re.search(r'(\d+)', hit)
    return int(m.group(1)) if m else 10**12


# clustered variants
# >GB00002469.1|6ccb96db9d81f6f110fe7cb5be5bdf5b|V_15083|S_241|P_11|r_0.143077314152|100.0
i = 0
vars_clusters = {}
for line in openfile(clustered_vars, 'r'):
    ch = line[0]
    if ch == '>':
        title = line[1:].strip()
        parts = title.split('|')
        md5_var = parts[1]
        vars_clusters[md5_var] = parts[0]
        i += 1

print("Processed sequence variants: " + str(i))

# raw fasta from samples
# >GB01020442S|An_2019_1acp_Bact|SRR5920425.4|POS=5|POS=253
# TACGTAGGGAGC...
vars_samples = {}
sample_occurence = {}
sample_id = 1
titleRead = False
i = 0
n = 0
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
            ###
            if vars_samples.has_key(seq):
                samples = vars_samples[seq]
                if samples.has_key(sample_name):
                    samples[sample_name] += 1
                else:
                    samples[sample_name] = 1
                    # count occurence
                    if sample_occurence.has_key(sample_name):
                        sample_occurence[sample_name] += 1
                    else:
                        sample_occurence[sample_name] = 1
                vars_samples[seq] = samples
            else:
                samples = {}
                samples[sample_name] = 1
                # count occurence
                if sample_occurence.has_key(sample_name):
                    sample_occurence[sample_name] += 1
                else:
                    sample_occurence[sample_name] = 1
                vars_samples[seq] = samples

print("Processed raw fasta from samples: " + str(i))

# sort sample names based on occurence
sample_names = create_sample_names(sample_occurence)

# SAVE PAIRING TABLE
i = 0
fp = open("VARIANTS_TABLE_SAMPLE_PAIRS.txt", 'w')
fp.write('SampleName\tSampleID\toccurence\n')
for sample_name in sample_names:
    fp.write(str(sample_names[sample_name]) + '\t' + sample_name + '\t' + str(sample_occurence[sample_name]) +'\n')
    i += 1
fp.close()

print("Sample pairs were saved: " + str(i))


# CREATE TABLE IF NOT EXISTS `variants` (
#   `hash` varchar(32) NOT NULL,
#   `samples` MEDIUMTEXT NOT NULL,
#   `abundances` MEDIUMTEXT NOT NULL,
#   `cluster` varchar(12) NOT NULL,
#   `sequence` TEXT NOT NULL
# );
i = 0
out_file = gzip.open("VARIANTS_TABLE.txt.gz", 'wt')
for seq in vars_samples:
    md5_var = hashlib.md5(seq.encode()).hexdigest()
    samples = vars_samples[seq]
    ids_str, counts_str = samples_to_strings(samples)
    if vars_clusters.has_key(md5_var):
        out_file.write(md5_var + '\t' + ids_str + '\t' + counts_str + '\t' + str(hit_numeric_key(vars_clusters[md5_var])) + '\t' + seq +'\n')
    else:
        out_file.write(md5_var + '\t' + ids_str + '\t' + counts_str + '\t' + '-' + '\t' + seq +'\n')
    i += 1
out_file.close()

print("Variants created: " + str(i))
print("DONE :)")
