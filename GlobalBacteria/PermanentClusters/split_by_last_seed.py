__author__ = 'vetrot'

import sys
import os
import gzip

input_fasta = sys.argv[1]       # GB_VOL1_PHYLUM_p__Latescibacterota.fas
seeds_text = sys.argv[2]        # GB_VOL1_PHYLUM_p__Latescibacterota.fas.97.0.seed_seqs

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################

fp = open(input_fasta + ".input_seeds_seqs", 'w')

lines = openfile(seeds_text).readlines()

# Write all except last line
for line in lines[:-1]:
    fp.write(line)

fp.close()

# store last line
last_seed = lines[-1].rstrip()

print("LAST SEED")
print(last_seed)
print("")


fp1 = open(input_fasta + ".done", 'w')
fp2 = open(input_fasta + ".undone", 'w')
done = True
filled = False
for n, line in enumerate(openfile(input_fasta)):
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line.rstrip()
            filled = True
    if filled:
        if seq == last_seed:
            done = False
        if done:
            fp1.write(">" + title + '\n')
            fp1.write(seq + '\n')
        else:
            fp2.write(">" + title + '\n')
            fp2.write(seq + '\n')
        filled = False
fp1.close()
fp2.close()
print("done")