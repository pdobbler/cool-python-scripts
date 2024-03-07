__author__ = 'vetrot'

import sys
import gzip
import hashlib

ITS_fasta = sys.argv[1]
ITS_clust_and_binned = sys.argv[2] # REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz.qualified.seeds
out_tab = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

clustered = {}
for line in openfile(ITS_clust_and_binned):
    ch = line[0]
    if ch == '>':
        title = line[1:].strip()
        hash_code =  line[1:].strip().split('|')[1]
        if clustered.has_key(hash_code):
            print("Error: multiple hash codes!!! " + hash_code)
        else:
            clustered[hash_code] = 0

print("Clustered variants loaded...")

samples_in = {}
samples_out = {}
for line in openfile(ITS_fasta):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            hash_code = hashlib.md5(line.strip().encode()).hexdigest()
            sample = title.split('|')[0]
            if not samples_in.has_key(sample):
                samples_in[sample] = 0
                samples_out[sample] = 0
            # test in/out
            if clustered.has_key(hash_code):
                samples_in[sample] += 1
            else:
                samples_out[sample] += 1

print("Sequnces processed...")

fp = open(out_tab, 'w')
fp.write('sample\tfound\tnot-found\n')
for sample in samples_in:
    fp.write(sample+'\t'+str(samples_in[sample])+'\t'+str(samples_out[sample])+'\n')
fp.close()

print("Done :]")

