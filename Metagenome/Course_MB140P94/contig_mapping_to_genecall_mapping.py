__author__ = 'vetrot'

import sys
import os

gene_call_fasta = sys.argv[1]
##title:
#>k141_20_1_453_+
ctg_mapping_tab = sys.argv[2]
##contig
#k141_43039
header = ''
abundances = {}
for n, line in enumerate(open(ctg_mapping_tab)):
    if n == 0:
        header = line.rstrip()
    else:
        vals = line.rstrip().split('\t')
        line_vals = ''
        for x in range(1,len(vals)):
            line_vals = line_vals + '\t'+vals[x]
        abundances[vals[0]] = line_vals

print("mapping table read...")

title = ''
sequence = ''
filled = False
genes_names = {}
for n, line in enumerate(open(gene_call_fasta)):
    if n % 20000 == 0:
        print(n / 2)
    if n % 2 == 0:
        title = line.rstrip()
        #print title
        if title[0] != '>':
            print("fasta format error...")
            break
    else:
        if n % 2 == 1:
            sequence = line.rstrip()
            filled = True

    if filled:
        tp = title[1:].rsplit('_',3)
        genes_names[title[1:]] = tp[0]
        filled = False

print("genecall fasta read...")

fp = open(ctg_mapping_tab + "_genecall.txt", 'w')
fp.write(header + "\n")
for name in genes_names:
    new_line = name + abundances[genes_names[name]]
    fp.write(new_line + "\n")
fp.close()

print("Done :)")

