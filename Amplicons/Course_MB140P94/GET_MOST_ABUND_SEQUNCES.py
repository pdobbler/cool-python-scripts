__author__ = 'vetrot'

import sys
import operator
import collections

fasta_file = sys.argv[1]
fasta_most = sys.argv[2]

print("Loading sequences...")

groups = {}
otu_sums = {}
filled = False
for n, line in enumerate(open(fasta_file)):
    if n % 200000 == 0:
        print(str(n / 2))
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            sequence = line.rstrip()
            filled = True
    if filled:
	parts = title.split('|')
        otu_name = parts[len(parts)-1]
        # sum otus...
        if otu_sums.has_key(otu_name):
            otu_sums[otu_name] = otu_sums[otu_name] + 1
        else:
            otu_sums[otu_name] = 1
        # attach sequences...
        if groups.has_key(otu_name):
            seqs = groups[otu_name]
            if seqs.has_key(sequence):
                seqs[sequence] = seqs[sequence] + 1
            else:
                seqs[sequence] = 1
            groups[otu_name] = seqs
        else:
            seqs = {sequence: 1}
            groups[otu_name] = seqs
        filled = False

# sort otus sums...
sorted_otu_sums = sorted(otu_sums.items(), key=operator.itemgetter(1), reverse=True)
print("Searching the most abundant sequences - please wait...")

# write fasta...
# CL0001|MOSTABUND|n=15/1
fp = open(fasta_most, "w")
for otu_name in sorted_otu_sums:
    seqs = groups[otu_name[0]]
    sorted_seqs = sorted(seqs.items(), key=operator.itemgetter(1), reverse=True)
    if otu_name[1]>1:
        fp.write(">" + otu_name[0] + "|MOSTABUND|n=" + str(otu_name[1]) + "/" + str(sorted_seqs[0][1]) + "\n")
    else:
        fp.write(">" + otu_name[0] + "|SINGLETON|n=" + str(otu_name[1]) + "\n")
    fp.write(sorted_seqs[0][0] + "\n")


fp.close()
print("Done :]")
