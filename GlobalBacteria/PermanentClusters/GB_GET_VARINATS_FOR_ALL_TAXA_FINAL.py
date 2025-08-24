__author__ = 'vetrot'

import sys
import gzip

variants_table = sys.argv[1] # TABLE_PROCESSED.txt
taxonomy_table = sys.argv[2] # TAXONOMY_utf8.txt
output_dir = sys.argv[3]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################
# TAXONOMY
#############################################

#CREATE TABLE IF NOT EXISTS `clusters_tax` (
#0  `id` int NOT NULL PRIMARY KEY,
#1  `cluster` VARCHAR(12) NOT NULL,
#2  `Species` varchar(64) NOT NULL,
#3  `Genus` varchar(32) NOT NULL,
#  `sim` FLOAT,
#  `cov` FLOAT,
#  `full_tax` TEXT NOT NULL,
#  `hash` varchar(32) NOT NULL
#);

# 26      GB00000026.1    -       -       100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__; s__     c0d5395792eadbf5f62e8ffb14fa0262
# 1       GB00000001.1    Bradyrhizobium sp000244915      Bradyrhizobium_503372   100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rhizobiales_505101; f__Xanthobacteraceae; g__Bradyrhizobium_503372; s__Bradyrhizobium sp000244915    6472eb8b1e09f892aca2f23182962903

sh_names = {}
sp_names = {}
gen_names = {}

for line in openfile(taxonomy_table):
    vals = line.rstrip().split('\t')
    sh_names[vals[0]] = vals[1]
    sp_names[vals[0]] = vals[2].replace(" ", "_")
    gen_names[vals[0]] = vals[3].replace(" ", "_")

sh_vars = {}
sp_vars = {}
gen_vars = {}

hash = {}
sequences = {}
markers = {}
samples = {}
abundances = {}
index = 0

# VARIANTS_TABLE.txt.gz
#CREATE TABLE IF NOT EXISTS `variants` (
#0  `hash` varchar(32) NOT NULL,
#1  `samples` MEDIUMTEXT NOT NULL,
#2  `abundances` MEDIUMTEXT NOT NULL,
#3  `BG` int NOT NULL,
#4  `sequence` TEXT NOT NULL
#);

# 6c8c16ea9ff8a0964d2444f30ddb5650        5404;1835;2061  1;1;2   218777  TACAGAGGGTGCGAGCGTTGTCCGGATTTATTGGGCGTAAAGAGCGTGTAGGCGGTTCGGTAGGTCCGTTGTGAAAACTCGAGGCTCAACCTCGAGACGCCGATGGAAACCCCCGAACTAGAGTCCGGAAGAGGAGAGTGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAAGAACACCCGTGGCTAAGGCGGCTCTCTAGTACGGTACTGACGCTGAGACGCGAAAGCGTGGGGAGCGAACAGG
# 262801bba46ead85dea520ae0234559b        6887    1       -       TACAGAGGTGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGCCGCCTAAGTCAGACGTGAAATCCCCCGGCTTAGCCTGGGAACTGCGTCTGATACTGGGTGGCTTGAGTTCGGGAGAGGGATGCGGAATTCCAGGTGTAGCGGTGAAATGCGTAGATATCTGGGGGAACACCGGTGGCGAAGGCGGCATCCTGGACCGAAACTGACGCTGAGGCGCGAAAGCTAGGGGAGCAAACGGG

for line in openfile(variants_table):
    vals = line.rstrip().split('\t')
    samples[index] = vals[1].split(';')
    abundances[index] = vals[2].split(';')
    if vals[3] != '-': # cluste is defined
        hash[index] = vals[0]
        sequences[index] = vals[4]
        #######################
        cl_name = sh_names[vals[3]]
        if sh_vars.has_key(cl_name):
            sh_vars[cl_name].append(index)
        else:
            sh_vars[cl_name] = []
            sh_vars[cl_name].append(index)

        sp_name = sp_names[vals[3]]
        if sp_name != '-':
            if sp_vars.has_key(sp_name):
                sp_vars[sp_name].append(index)
            else:
                sp_vars[sp_name] = []
                sp_vars[sp_name].append(index)

        gen_name = gen_names[vals[3]]
        if gen_name != '-':
            if gen_vars.has_key(gen_name):
                gen_vars[gen_name].append(index)
            else:
                gen_vars[gen_name] = []
                gen_vars[gen_name].append(index)
        #######################
    index += 1

# compute sample ITS sums...
total = {}
for index in samples:
    for i in range(len(samples[index])):
        s = samples[index][i]
        a = abundances[index][i]
        if total.has_key(s):
            total[s] = total[s] + int(a)
        else:
            total[s] = int(a)
#for s in total:
#    print("Sample "+s+" has "+ str(total[s]) +" sequnces.")
print("sample total sums were computed...")

# generate files...
print("output: " + output_dir + " #of annotated variants is " + str(index))

print(" #genera " + str(len(gen_vars)))
print(" #species " + str(len(sp_vars)))
print(" #SH " + str(len(sh_vars)))
# sh
for name in sh_vars:
    i = 1
    fp = gzip.open(output_dir + "SH_" + name+".fas", "wb")
    for index in sh_vars[id]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash[index] + "|SampleID_" + s + "|sh_" + name + "|marker_16S|abund_"  + a +"_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()

print("SH variants were written...")

# sp
for name in sp_vars:
    if "_sp." not in name:
        i = 1
        fp = gzip.open(output_dir + "species_" + name+".fas", "wb")
        for index in sp_vars[id]:
            for i in range(len(samples[index])):
                s = samples[index][i]
                a = abundances[index][i]
                fp.write(">" + hash[index] + "|SampleID_" + s + "|species_" + name + "|marker_16S|abund_" + a + "_total_" + str(total[s]) + "\n")
                fp.write(sequences[index] + "\n")
                i += 1
        fp.close()
    else:
        print("Ignored: "+name)

print("Species variants were written...")

# gen
for name in gen_vars:
    i = 1
    fp = gzip.open(output_dir + "genus_" + name+".fas", "wb")
    for index in gen_vars[id]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash[index] + "|SampleID_" + s + "|genus_" + name + "|marker_16S|abund_" + a + "_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()

print("Genus variants were written...")

########
print("Done :)")


