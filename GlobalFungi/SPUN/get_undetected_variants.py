__author__ = 'vetrot'

import sys
import os
import gzip

new_vars_fasta_file = sys.argv[1]
old_fasta_file = sys.argv[2]
unfind_vars_out = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

variants = {}
filled = False
for n, line in enumerate(openfile(new_vars_fasta_file)):
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line.rstrip()
            filled = True
    if filled:
        variants[seq] = title
        filled = False
print("New variants loaded... "+str(len(variants)))

filled = False
for n, line in enumerate(openfile(old_fasta_file)):
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line.rstrip()
            filled = True
    if filled:
        if variants.has_key(seq):
            variants[seq] = ""
        filled = False
print("Sequences checked... ")

fp = open(unfind_vars_out, 'w')
for seq in variants:
    if not variants[seq] == "":        
        fp.write(">" + variants[seq] + '\n')
        fp.write(seq + '\n')
fp.close()
print("Done :]")

