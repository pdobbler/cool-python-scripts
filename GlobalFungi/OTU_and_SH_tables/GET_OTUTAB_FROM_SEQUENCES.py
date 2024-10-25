__author__ = 'vetrot'

import sys
import operator
import collections

fasta_file = sys.argv[1]
out_tab = sys.argv[2]


#>GF05023557S|Fernan_2020_JZ12|e7c017be559f77ec|OTU0027655
# fill table...
tab_dict = {}
samples_sums = {}
otu_sums = {}
filled = False
for line in open(fasta_file):
    if line[0] == ">":
        vals = line.rstrip().split('|')
        # sum otus...
        otu_name = vals[3]
        sample_name = vals[0]

        if otu_sums.has_key(otu_name):
            otu_sums[otu_name] = otu_sums[otu_name] + 1
        else:
            otu_sums[otu_name] = 1
        # sum samples...
        if samples_sums.has_key(sample_name):
            samples_sums[sample_name] = samples_sums[sample_name] + 1
        else:
            samples_sums[sample_name] = 1
        # fill dict...
        if tab_dict.has_key(otu_name):
            samples_counts = tab_dict[otu_name]
            if samples_counts.has_key(sample_name):
                samples_counts[sample_name] = samples_counts[sample_name] + 1
            else:
                samples_counts[sample_name] = 1
            tab_dict[otu_name] = samples_counts
        else:
            samples_counts = {sample_name: 1}
            tab_dict[otu_name] = samples_counts

# sort otus sums...
sorted_otu_sums = sorted(otu_sums.items(), key=operator.itemgetter(1), reverse=True)
# sort samples sums...
sorted_samples_sums = collections.OrderedDict(sorted(samples_sums.items()))

# write table...
fp = open(out_tab, "w")
# header
line = "OTU"
for sample in sorted_samples_sums:
    print("Sample " + sample + " has " + str(sorted_samples_sums[sample]) + " sequences")
    line += "\t" + sample
fp.write(line + "\n")
# counts
for otu in sorted_otu_sums:
    line = otu[0]
    #print("Otu " + str(otu[0]))
    for sample in sorted_samples_sums:
        samples_counts = tab_dict[otu[0]]
        val = 0
        if samples_counts.has_key(sample):
            val = samples_counts[sample]
        line += "\t" + str(val)
    fp.write(line + "\n")

fp.close()
print("Done :]")



