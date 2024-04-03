__author__ = 'vetrot'

import sys
import operator
import os

fasta_file = sys.argv[1]
out_fasta = sys.argv[2]

#################
# DEREPLICATION #
#################

cmd = 'usearch -fastx_uniques '+fasta_file+' -fastaout '+fasta_file+'.uniques.fa -sizeout -relabel Uniq'
print(cmd)

stream = os.popen(cmd)
output = stream.read()
print(output)

print("Uniques sequences were generated using fastx_uniques command...")

##############
# CLUSTERING #
##############

cmd = 'usearch -cluster_otus '+fasta_file+'.uniques.fa -minsize 1 -otus '+fasta_file+'.otus.fa -relabel Otu -uparseout '+fasta_file+'.uparse.txt'
print(cmd)

stream = os.popen(cmd)
output = stream.read()
print(output)

print("Clustering on 97 similarity level was performed using cluster_otus command...")

####################
# LABELING UNIQUES #
####################
cmd1 = r'/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}';
cmd = "awk '"+cmd1+"' "+fasta_file+".uniques.fa > "+fasta_file+".uniques.linear.fa"
print(cmd)

stream = os.popen(cmd)
output = stream.read()
print(output)

print("Robert Edgar's linebreaks in "+fasta_file+".uniques.fa were removed by awk...")

unique_seqs = {}
filled = False
for n, line1 in enumerate(open(fasta_file+".uniques.linear.fa")):
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
for line in open(fasta_file+'.uparse.txt'):
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

print("OTU labels were extracted from uparseout file...")

##################
# LABELING INPUT #
##################

fp = open(out_fasta, "w")
fpc = open(out_fasta + ".chimeric", "w")
filled = False
for n, line1 in enumerate(open(fasta_file)):
    if n % 200000 == 0:
        print(str(n / 2))
    if n % 2 == 0:
        header = line1.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line1.rstrip()
            filled = True
    if filled:
        uniq_title = unique_seqs[seq]
        label = dict[uniq_title]
        if label[:3] != 'otu':
            fpc.write(">" + header + "|" + label + "\n")
            fpc.write(seq + "\n")
        else:
            #change otu number to fixed digit text...
            otu_name = "OTU"
            for i in range(max_len - len(label)):
                otu_name += "0"
            otu_name += label[3:]
            #write sequence...
            fp.write(">" + header + "|" + otu_name + "\n")
            fp.write(seq + "\n")
        filled = False
fp.close()
fpc.close()

print("Sequences were renamed based on the generated OTUs....")

print("All is done!")

