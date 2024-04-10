__author__ = 'vetrot'

import sys
import os

annotation = sys.argv[1]
short_to_tax = sys.argv[2]
output_file = sys.argv[3]
outtype = sys.argv[4]           # D/

if outtype == 'D':
    print("Out-type is: DEFAULT")
if outtype == 'U':
    print("Out-type is: UPDATE")

#fwd_k141_10000352_1_282_-       jgi|Dacma1|778102|KOG1368|4.1.2.5       52.5    80      38      0       13      92      20      99      5.1e-15 87.0

names = {}
for line in open(short_to_tax):
    vals = line.rstrip().split('\t')
    names[vals[0]] = vals[1]

print('Taxon paires were loaded - '+str(len(names)))

fp = open(output_file, "w")
for line in open(annotation):
    vals = line.rstrip().split('\t')
    new_line = vals[0]
    for i in range(1,len(vals)):
        if i == 1:
            tax = vals[i].split('|')[1]
            if names.has_key('['+tax+']'):
                new_line = new_line + '\t' + names['['+tax+']']
            else:
                print('Taxa abreciation was not found - '+tax)
                new_line = new_line + '\t' + vals[i]
            if outtype == U:
                new_line += '|' + vals[i].split('|')[3]
        else:
            new_line = new_line + '\t' + vals[i]
    #print(new_line)
    fp.write(new_line + '\n')

fp.close()

print('Done :]')



