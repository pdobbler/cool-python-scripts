__author__ = 'vetrot'

import sys

otu_tab = sys.argv[1]
fasta_file = sys.argv[2]

fp = open(fasta_file, 'w')
sample_names = {}
index = 0
for n, line in enumerate(open(otu_tab)):
    vals = line.rstrip().split('\t')
    if n<1:
        i=1
        while i < len(vals):
            sample_names[i] = vals[i]
            #sample_sums[i] = 0.0
            print vals[i]
            i = i + 1
    else:
        i=1
        while i < len(vals):
            #sample_sums[i] = sample_sums[i] + float(vals[i])/float(copy_numbers[vals[0]])
            #print str(sample_sums[sample_names[i-1]])
            n = 0
            while n < int(vals[i]):
                index = index +1
                line_new = str(index) + "|" + sample_names[i] + "|" + vals[0]
                #print line_new
                fp.write(">"+line_new + "\n")
                fp.write("ACGT" + "\n")
                n = n + 1
            i = i + 1
fp.close()

print "Done"