__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

input_file = sys.argv[1]                # OTU_TAB_SOIL_NOSINGLE_FUNGI.txt
output_file = sys.argv[2]    # OTU_TAB_SOIL_NOSINGLE_FUNGI_TRANSPOSED.txt

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# Transpose a tab-separated file
#############################################

# Read the input file
with open(input_file, 'r') as infile:
    data = [line.strip().split('\t') for line in infile]

# Transpose the data
transposed_data = list(map(list, zip(*data)))

# Write the transposed data to the output file
with open(output_file, 'w') as outfile:
    for row in transposed_data:
        outfile.write('\t'.join(row) + '\n')

print("DONE :]")
