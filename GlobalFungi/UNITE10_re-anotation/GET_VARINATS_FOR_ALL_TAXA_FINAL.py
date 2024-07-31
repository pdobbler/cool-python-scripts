__author__ = 'vetrot'

import sys

variants_table = sys.argv[1] # TABLE_PROCESSED.txt
taxonomy_table = sys.argv[2] # TAXONOMY_utf8.txt
output_dir = sys.argv[3]

#CREATE TABLE IF NOT EXISTS `taxonomy` (
#0  `SH` varchar(32) NOT NULL,
#1  `Ecology` varchar(64) NOT NULL,
#2  `Kingdom` varchar(64) NOT NULL,
#3  `Phylum` varchar(64) NOT NULL,
#4  `Class` varchar(64) NOT NULL,
#5  `Order` varchar(64) NOT NULL,
#6  `Family` varchar(64) NOT NULL,
#7  `Genus` varchar(64) NOT NULL,
#8  `Species` varchar(64) NOT NULL,
#9  `genus_id` int NOT NULL,
#10  `species_id` int NOT NULL,
#11 `SH_id` int NOT NULL
#);

sh_names = {}
sp_names = {}
gen_names = {}

for line in open(taxonomy_table):
    vals = line.rstrip().split('\t')
    sh_names[vals[11]] = vals[0]
    sp_names[vals[10]] = vals[8].replace(" ", "_")
    gen_names[vals[9]] = vals[7].replace(" ", "_")

sh_vars = {}
sp_vars = {}
gen_vars = {}

hash = {}
sequences = {}
markers = {}
samples = {}
abundances = {}
index = 0

for line in open(variants_table):
    vals = line.rstrip().split('\t')
    samples[index] = vals[1].split(';')
    abundances[index] = vals[2].split(';')
    if vals[4] != '0':
        hash[index] = vals[0]
        sequences[index] = vals[7]
        markers[index] = vals[3]
        #######################
        if sh_vars.has_key(vals[4]):
            sh_vars[vals[4]].append(index)
        else:
            sh_vars[vals[4]] = []
            sh_vars[vals[4]].append(index)

        if vals[5] != '0':
            if sp_vars.has_key(vals[5]):
                sp_vars[vals[5]].append(index)
            else:
                sp_vars[vals[5]] = []
                sp_vars[vals[5]].append(index)

        if vals[6] != '0':
            if gen_vars.has_key(vals[6]):
                gen_vars[vals[6]].append(index)
            else:
                gen_vars[vals[6]] = []
                gen_vars[vals[6]].append(index)
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
for id in sh_vars:
    name = sh_names[id]
    i = 1
    fp = open(output_dir + "SH_" + name+".fas", "w")
    for index in sh_vars[id]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash[index] + "|SampleID_" + s + "|sh_" + name + "|marker_" + markers[index] + "|abund_"  + a +"_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()

print("SH variants were written...")

# sp
for id in sp_vars:
    name = sp_names[id]
    if "_sp." not in name:
        i = 1
        fp = open(output_dir + "species_" + name+".fas", "w")
        for index in sp_vars[id]:
            for i in range(len(samples[index])):
                s = samples[index][i]
                a = abundances[index][i]
                fp.write(">" + hash[index] + "|SampleID_" + s + "|species_" + name + "|marker_" + markers[index] + "|abund_" + a + "_total_" + str(total[s]) + "\n")
                fp.write(sequences[index] + "\n")
                i += 1
        fp.close()
    else:
        print("Ignored: "+name)

print("Species variants were written...")

# gen
for id in gen_vars:
    name = gen_names[id]
    i = 1
    fp = open(output_dir + "genus_" + name+".fas", "w")
    for index in gen_vars[id]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash[index] + "|SampleID_" + s + "|genus_" + name + "|marker_" + markers[index] + "|abund_" + a + "_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()

print("Genus variants were written...")

########
print("Done :)")


#CREATE TABLE IF NOT EXISTS `variants` (
#0  `hash` varchar(32) NOT NULL,
#1  `samples` MEDIUMTEXT NOT NULL,
#2  `abundances` MEDIUMTEXT NOT NULL,
#3  `marker` varchar(4) NOT NULL,
#4  `SH` int NOT NULL,
#5  `species` int NOT NULL,
#6  `genus` int NOT NULL,
#7  `sequence` TEXT NOT NULL
#);

