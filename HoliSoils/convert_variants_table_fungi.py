__author__ = 'vetrot'

import sys
import gzip

variants_file     = sys.argv[1]  # GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz
samples_to_filter = sys.argv[2]  # GF5_samples_holisoils.txt

def openfile(filename, mode='r'):
    # For Py3 use 'rt'/'wt' with encoding if needed.
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# -----------------------------
# Load Holisoils samples (map ID -> NAME)
# File format example: "GF01000413S\t574"
# -----------------------------
samples = {}         # key: sample_id (string), value: sample_name
samples_found = {}   # key: sample_name, value: 0/1
for line in openfile(samples_to_filter, 'r'):
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    vals = line.split('\t')
    if len(vals) < 2:
        continue
    sname, sid = vals[0], vals[1]
    samples[sid] = sname
    samples_found[sname] = 0
print("Samples to filter loaded: " + str(len(samples)))

# -----------------------------
# Prepare outputs
# NOTE: This writes 5 fields to variants:
#   var_id, cl_id, hash, sequence, marker
# -----------------------------

# FUNGAL
# CREATE TABLE IF NOT EXISTS `variants` (
#   `hash` varchar(32) NOT NULL,
#   `samples` MEDIUMTEXT NOT NULL,
#   `abundances` MEDIUMTEXT NOT NULL,
#   `marker` varchar(4) NOT NULL,
#   `SH` int NOT NULL,
#   `sequence` TEXT NOT NULL
# );
# 8287dcc6431a451c2386c7f7d0131c04        14      1       ITS2    0       CAACCATCAAGCCCTTGCTTGTGTTGGGGTCCTGCGGCTGCCCGCAGGCCCTGAAAACCAGTGGCGGGCTCCCGAGTCACACCGAGGGCAGTAATACATCTCGCTTTGGTCGGGGGGGGGGGTCCTGCCGGTAAAAACCCCCCTTTCCAAA
# 90e35b1036fc89326030639387d2fc4a        15      1       ITS1    0       TTGATTTTTAAAGATGTGCTGGCATTAATGCATGTGCACTCTTCACAAAACCAATATCCACCTGTGCACACTCTGTAAGCAACAGGGTCTTTACCCTTTGCTTATGTCCCTTTACAACCAACCATTGATTTATTTGACTGTAATTTTGAAATAACTAAAA
# 660f8b39108cc20d0dbacee2e154f980        16;17   3;3     ITS2    4       ACACCTCAACTCTTCATGGTTTTCCATGATGAGCTTGGACTTTGGGGGTCTTGCTGGCCTGCGGTCGGCTCCTCTCAAATGAATCAGCTTGCCAGTGTTTGGTGGGCATCACAGGTGTGATAACTATCTACGCTTGTGGTTTTCCACCAGGTAACCTTCAGCAGTGGAGGTTCGCTGGAGCTCACAGATGTCTCTCCTCAGTGAGGGCAGCCCTTTGTAT


# CREATE TABLE IF NOT EXISTS `variants` (
#   `id` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL,
#   `hash` varchar(32) NOT NULL,
#   `sequence` TEXT NOT NULL
#   `marker` varchar(4) NOT NULL,
# );

# CREATE TABLE IF NOT EXISTS `samplevar` (
#   `id` bigint(20) unsigned NOT NULL,
#   `variant` int(10) unsigned NOT NULL,
#   `sample` int(10) unsigned NOT NULL,
#   `abundance` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL
# );


fp_variants  = gzip.open("VARIANTS_FUN_variants.txt.gz",  'wb')
fp_samplevar = gzip.open("VARIANTS_FUN_samplevar.txt.gz", 'wb')

sv_id  = 1
var_id = 1

for line in openfile(variants_file, 'r'):
    line = line.rstrip('\n')
    if not line or line.startswith('#'):
        continue

    vals = line.split('\t')
    if len(vals) < 6:
        # Unexpected row, skip
        continue

    md5 = vals[0]
    sn  = vals[1].split(';')  # sample IDs (strings)
    sa  = vals[2].split(';')  # abundances (strings)
    marker = vals[3]
    cl_id  = vals[4]
    seq    = vals[5]

    # Guard length match
    if len(sn) != len(sa):
        # Skip malformed line
        continue

    # Indices of samples that are in the Holisoils list
    keep_idx = [i for i, sid in enumerate(sn) if sid in samples]

    if keep_idx:
        # Write one variants row
        fp_variants.write(("{0}\t{1}\t{2}\t{3}\t{4}\n")
                          .format(var_id, cl_id, md5, seq, marker))

        # Write only the kept samplevar rows
        for i in keep_idx:
            sid = sn[i]          # sample ID
            abund = sa[i]
            sname = samples[sid] # map ID -> name for found tracking
            samples_found[sname] = 1

            fp_samplevar.write(("{0}\t{1}\t{2}\t{3}\t{4}\n")
                               .format(sv_id, var_id, sid, abund, cl_id))
            sv_id += 1

        var_id += 1

fp_samplevar.close()
fp_variants.close()

print("Tables saved...")

# Report missing samples by NAME (consistent with init)
for sname, seen in samples_found.items():
    if seen == 0:
        print("WARNING: sample " + sname + " was not found!")

print("Done :]")
