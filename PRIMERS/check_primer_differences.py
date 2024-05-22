__author__ = 'vetrot'

import sys
import datetime
from Bio import Align

original_primers = sys.argv[1]      # could be empty
new_primers = sys.argv[2]
minMM = int(sys.argv[3])
output_primers = sys.argv[4]

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
#aligner.match_score = 1
#aligner.mismatch_score = -1

print("Pairwise Aligner mode: " + aligner.mode)
print("Pairwise Aligner open_gap_score: " + str(aligner.open_gap_score))
print("Pairwise Aligner extend_gap_score: " + str(aligner.extend_gap_score))
print("Pairwise Aligner match_score: " + str(aligner.match_score))
print("Pairwise Aligner mismatch_score: " + str(aligner.mismatch_score))

def compare_strings(str1, str2):
    alignments = aligner.align(str1, str2)
    alignment_str = str(alignments[0])
    match_line = alignment_str.split('\n')[1]
    mm = len(match_line)-match_line.count('|')
    return mm


PRIMER_LIB = {}

# 528F_H001 TGTATAGCGGTAATTCCAGCTCCAA
PASSED_LIST = []
for line in open(original_primers):
    vals = line.rstrip().split('\t')
    if vals[0] != "":
        if vals[0] in PRIMER_LIB:
            print("ERROR: DUPLICATE PRIMER NAME - " + vals[0])
            sys.exit()
        else:
            PRIMER_LIB[vals[0]] = vals[1]
            PASSED_LIST.append(vals[0])

# 528F_7    TAATACTAGCGGTAATTCCAGCTCCAA
NEW_LIST = []
for line in open(new_primers):
    vals = line.rstrip().split('\t')
    if vals[0] != "":
        if vals[0] in PRIMER_LIB:
            print("ERROR: DUPLICATE PRIMER NAME - " + vals[0])
            sys.exit()
        else:
            PRIMER_LIB[vals[0]] = vals[1]
            NEW_LIST.append(vals[0])

GOOD_LIST = {}
i = 0
while NEW_LIST:
    top_item = NEW_LIST.pop(0)
    top_seq = PRIMER_LIB[top_item]
    i += 1
    print(":::ROUND: "+ str(i))
    skip = False
    for name in PASSED_LIST:
        test_seq = PRIMER_LIB[name]
        print("****************TESTING******")
        print(test_seq + "   " +name)
        print(top_seq + "   " +top_item)
        mm = compare_strings(test_seq, top_seq)
        print(mm)
        if mm < minMM:
            skip = True
            print("Too similar skipping...")
            break
    if not skip:
        print("GOOD SEQ>" + top_item)
        GOOD_LIST[top_item] = top_seq
        PASSED_LIST.append(top_item)


# save the good ones
out_file = open(output_primers, "w")
for name in GOOD_LIST:
    out_file.write(name + "\t" + GOOD_LIST[name] + "\n")
out_file.close()

print("Done :)")




