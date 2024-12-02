#!/usr/bin/env python3

__author__ = 'thiagp'

import os
import sys
from collections import defaultdict
import argparse

#in_cazy = sys.argv[1]
out_cazy = sys.argv[2]

#parser = argparse.ArgumentParser()

#parser.add_argument("-i", "in_table", help="all_CAZy table with all res. from hmm.out cated")
#parser.add_argument("-o", "out_table", help="output file name with one gene per line with multiple CAZy")
#parser.add_argument("-l", "--label", help="label to be added to output file", type=str, default = '')
#args=parser.parse_args()


gene_dict = defaultdict(list)
with open(sys.argv[1]) as in_cazy:
	for line in in_cazy:
		if line.strip() and not 'HMM Profile' in line:
			vals = line.rstrip().split('\t')
			#print(vals)
			gene = vals[2]
			cazy_name = vals[0].split('.', 1)[0]
			evals = vals[4]
			#(gene, *cazy) = (vals[3], vals[0].split('.', 1)[0])
			#(evalu, start) = (vals[5], vals[7])
			#(*codes, study) = line.split('\t')
			if gene in gene_dict:
				cazydict = gene_dict[gene]
				#print(cazydict)
				cazydict[cazy_name] = float(evals)
				#gene_dict[gene][cazy] = evalu, start
			else:
				cazydict = defaultdict(list)
				cazydict[cazy_name] = float(evals)
				gene_dict[gene] = cazydict
				#dictionary[gene][cazy] = evalu, start
	#print "raw all_CAZy table read"
	#

fp = open(out_cazy, 'w')
for gene in gene_dict:
	cazydict = gene_dict[gene]
	if len(cazydict) > 0:
		sort_orders = sorted(cazydict.items(), key=lambda x: x[1])
		cazylist =[]
		for i in cazydict.keys():
			cazylist.append(i)
		fp.write(gene + '\t' + str(sort_orders[0][1]) + '\t' +  ';'.join(cazylist) + '\n')
fp.close()
#
print("done?")


