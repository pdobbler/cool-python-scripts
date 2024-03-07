__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

ITS2_fasta = sys.argv[1] # PERMANENT_CLUSTERS/GF5_ALL_SAMPLES.fa.its2.gz
cluster_names = sys.argv[2] # TOP_2000_CLUSTER_NAMES.txt
clustered_variants = sys.argv[3] # FINAL/GF5_ALL_SAMPLES_ITS2_CLUSTANDBINNED.fa
out_tab = sys.argv[4]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
	if filename.endswith('.gz'):
		return gzip.open(filename, mode)
	else:
		return open(filename, mode)


#############################################
# EXTRACT ABUNDANCES
#############################################
clusters = {}
for line in openfile(cluster_names):
	vals = line.rstrip().split('\t')
	clusters[vals[0]] = {}

print("Cluster names loaded... (" + str(len(clusters)) + ")")

# >CL000001|e399dc7c0823811b06f9eea10fa4ea6c|V_15570922|S_31816|P_429|r_0.77938871414|100.0
# >CL192043|8645414f339ccc92e0cd08b6f081f469|V_1|S_1|P_1|r_0.00024882711719|BINNED_198.5
variants = {}
for line in openfile(clustered_variants):
        ch = line[0]
	if ch == '>':
		title = line[1:].strip()
		vals = title.split('|')
		if clusters.has_key(vals[0]):
			variants[vals[1]] = vals[0]

print("Variants were selected... (" + str(len(variants)) + ")")

# >GF05023293S|Lee123_2020_OO12|735e97b0c80a537f
# AAAAAAAAAAAAATATCGTTACATCTTTTTGGTGTTACGGATCTGGGTTATCCGGTTTTTAAGTCGGTTACCTAAAATTAAGGATTATATATAATGTGATACGTACTAAGATAAAGTCGTTAATCATTAAATTTATTACGTTATCAACTAATATGTTTAGTAGGTAGTGTATAATTTCTGA
i = 0
titleRead = False
sample_set = set()
for line in openfile(ITS2_fasta):
	ch = line[0]
	if ch == '>':
		titleRead=True
		title = line[1:].strip()
	else:
		if titleRead:
			titleRead=False
			i = i + 1
			variant = hashlib.md5(line.strip().encode()).hexdigest()
			if variants.has_key(variant):
				samples = clusters[variants[variant]]
				sample = title.split('|')[0]
				sample_set.add(sample)
				if samples.has_key(sample):
					samples[sample] += 1 
				else:
					samples[sample] = 1

print("Sequences were examined... (" + str(i) + ")")


sample_list = []
for sample in sample_set:
	sample_list.append(sample)

print("Sample set converted to list...")

# save it
fp = open(out_tab, "w")
newline = "CLUSTER"
for sample_name in sample_list:
    newline += "\t" + sample_name  
fp.write(newline + "\n")

for otu_name in clusters:
	newline = otu_name
	samples = clusters[otu_name]
	for sample_name in sample_list:
		if samples.has_key(sample_name):
			newline += "\t" + str(samples[sample_name])
		else:
			newline += "\t0"
	fp.write(newline + "\n")

fp.close()

print("DONE :]")

