__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

ITS2_fasta = sys.argv[1] # PERMANENT_CLUSTERS/GF5_ALL_SAMPLES.fa.its2.gz
clustered_variants = sys.argv[2] # FINAL/GF5_ALL_SAMPLES_ITS2_CLUSTANDBINNED.fa
out_tab = sys.argv[3]

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
# >CL000001|e399dc7c0823811b06f9eea10fa4ea6c|V_15570922|S_31816|P_429|r_0.77938871414|100.0
# >CL192043|8645414f339ccc92e0cd08b6f081f469|V_1|S_1|P_1|r_0.00024882711719|BINNED_198.5
# >CL000001|82780436adb026d0faf0ab5680b15b13|V_10502756|S_22243|P_404|r_358.89364931|100.0
variants = {}
for line in openfile(clustered_variants):
        ch = line[0]
	if ch == '>':
		title = line[1:].strip()
		vals = title.split('|')
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


# write table...
fp = open(clustered_variants + ".samples_to_otus_REDUCED.txt", "w")
fp.write("OTU\tsample\tabundance\n")
for otu_name in clusters:
	samples = clusters[otu_name]
	key_list = []
	value_list = []
	for key, value in samples.iteritems():
		key_list.append(key)
		value_list.append(value)
	fp.write(otu_name + "\t" + ";".join(key_list) + "\t" + ";".join(value_list) + "\n")
fp.close()

print("DONE :]")

