__author__ = 'vetrot'

import sys
from Bio import Align

text_file = sys.argv[1]
out_file = sys.argv[2]

aligner = Align.PairwiseAligner()
aligner.mode = "global"

print("Pairwise Aligner mode: " + aligner.mode)

fp = open(out_file, 'w')
for line in open(text_file):
	parts = line.strip().split('\t')
	alignments = aligner.align(parts[0], parts[1])
	score = alignments[0].score
	dl = len(parts[1]) - len(parts[0])
	fp.write(line.strip() + "\t" + str(score-len(parts[0])) + "\t"  + str(score-len(parts[1])) + "\t" + str(dl) + "\n")
fp.close()

print("DONE")


