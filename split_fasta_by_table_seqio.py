#!/usr/bin/env python3

__author__ = 'thiagp'

from Bio import SeqIO
import os
import sys
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("table", help="table with two columns, ID(same as seq header) and label to be split")
parser.add_argument("fasta", help="run fasta file to be split")
parser.add_argument("-l", "--label", help="label to be adde to output file", type=str, default = '')
args=parser.parse_args()


dictionary = defaultdict(list)

with open (sys.argv[1]) as studymap:
	for line in studymap.read().splitlines():
		(*codes, study) = line.split('\t')
		if study in dictionary:
			dictionary[study] += codes
		else:
			dictionary[study] = codes

print (dictionary) 


for i in dictionary:
	with open(i + '_' + args.label + '.fasta', 'w') as out:
		for seq_record in SeqIO.parse(sys.argv[2], "fasta"):
			samples = seq_record.id.split('|')[1]
			if samples in dictionary[i]:
				SeqIO.write(seq_record, out , "fasta-2line")
				



