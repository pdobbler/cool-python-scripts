__author__ = 'vetrot'

import sys
import operator
import os
import gzip

fasta_file = sys.argv[1] # GF5_ALL_SAMPLES.fa.its2.gz
out_fasta = sys.argv[2]
linear_fasta = sys.argv[3] # GF5_ALL_SAMPLES_ITS2_VARS_minsize2.fa
# >0e5ebdb4cd426f308dbf1307e041621f;size=27249999
# ATAGAGA...
# 
uparse_out = sys.argv[4]
# cf68cb2e631dd3a6c2f730720fec138e;size=13415448  otu10   dqt=66;
# 2ad7e1f00dec9a8dcb06330c5a38abbd;size=6329491   match   dqt=1;top=Otu14(99.6%);
# 5a2728462482027986193a018b81e8ee;size=48        perfect top=Otu47376(100.0%);

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# PROCESS
#############################################

unique_seqs = {}
filled = False
for n, line1 in enumerate(openfile(linear_fasta)):
    if n % 200000 == 0:
        print(str(n / 2))
    if n % 2 == 0:
        title = line1.rstrip()[1:]
    else:
        if n % 2 == 1:
            sequence = line1.rstrip()
            filled = True
    if filled:
        unique_seqs[sequence] = title
        filled = False

print("uniques were loaded to dictionary...")

dict = {}
max_len = 0
for line in openfile(uparse_out):
    vals = line.rstrip().split('\t')
    label = vals[1]
    if vals[1] == "match":
        label = vals[2].split(';')[1].split('=')[1].split('(')[0].lower()
        #print(vals[1]+label)
    if vals[1][:3] == "otu":
        label = vals[1]
        if len(vals[1]) > max_len:
            max_len = len(vals[1])
        #print(vals[1])
    if vals[1] == "perfect":
        label = vals[2].split('=')[1].split('(')[0].lower()
    # store it...
    if dict.has_key(vals[0]):
        print("Error - uniques are not uniques!!!")
    else:
        dict[vals[0]] = label

print("OTU labels were extracted from crazy uparseout file...")

##################
# LABELING INPUT #
##################

ch = 0
si = 0
cl = 0
fp = gzip.open(out_fasta, "wt")
fps = gzip.open(out_fasta + ".singletons", "wt")
fpc = gzip.open(out_fasta + ".chimeric", "wt")
filled = False
for n, line1 in enumerate(openfile(fasta_file)):
    if n % 200000 == 0:
        print(str(n / 2))
    if n % 2 == 0:
        header = line1.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line1.rstrip()
            filled = True
    if filled:
        if unique_seqs.has_key(seq):
            uniq_title = unique_seqs[seq]
            label = dict[uniq_title]
            if label[:3] != 'otu':
                fpc.write(">" + header + "|" + label + "\n")
                fpc.write(seq + "\n")
                ch += 1
            else:
                #change otu number to fixed digit text...
                otu_name = "OTU"
                for i in range(max_len - len(label)):
                    otu_name += "0"
                otu_name += label[3:]
                #write sequence...
                fp.write(">" + header + "|" + otu_name + "\n")
                fp.write(seq + "\n")
                cl += 1
        else:
            fps.write(">" + header + "\n")
            fps.write(seq + "\n")
            si += 1
        filled = False
fp.close()
fps.close()
fpc.close()

print("Sequences were renamed based on the generated OTUs - total: " + str(cl + si + ch))
print("Clustered to OTUs     : " + str(cl))
print("Removed chimeras      : " + str(ch))
print("Unclustered singletons: " + str(si))
