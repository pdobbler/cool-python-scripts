__author__ = 'vetrot'

import sys
import os
import gzip
import zipfile

variants_table = sys.argv[1] # GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz
taxonomy_table = sys.argv[2] # GF5_UNITE10_TAXONOMY_TABLE.txt.gz
output_dir = sys.argv[3]
samples_map = sys.argv[4] # GF5_RAW_TABLE_SAMPLES.txt.gz

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################

def zip_file(file_path):
    # Remove the .fas extension and add .zip extension
    zip_filename = os.path.splitext(file_path)[0] + '.zip'
    with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED, allowZip64=True) as zipf:
        # Add file to the zip file without directory paths
        zipf.write(file_path, os.path.basename(file_path))
    # Remove the original file
    os.remove(file_path)

#############################################

samples_mapping = {}
for line in openfile(samples_map):
    vals = line.rstrip().split('\t')
    samples_mapping[vals[0]] = vals[1]

print("Samples mapping loaded... (" + str(len(samples_mapping)) + ")")

#############################################

sh_names = {}

for line in openfile(taxonomy_table):
    vals = line.rstrip().split('\t')
    sh_names[vals[8]] = vals[0]+";"+vals[7].replace(" ", "_")+";"+vals[6].replace(" ", "_")

sh_vars = {}
sp_vars = {}
gen_vars = {}

hash_md5 = {}
sequences = {}
markers = {}
samples = {}
abundances = {}
index = 0

for line in openfile(variants_table):
    vals = line.rstrip().split('\t')
    samples[index] = vals[1].split(';')
    abundances[index] = vals[2].split(';')
    if vals[4] != '0':
        hash_md5[index] = vals[0]
        sequences[index] = vals[5]
        markers[index] = vals[3]
        # get taxa levels
        tax_levels = sh_names[vals[4]].split(';')
        sh_name = tax_levels[0]
        sp_name = tax_levels[1]
        gen_name = tax_levels[2]
        #sh
        if sh_vars.has_key(sh_name):
            sh_vars[sh_name].append(index)
        else:
            sh_vars[sh_name] = []
            sh_vars[sh_name].append(index)
        # species
        if sp_vars.has_key(sp_name):
            sp_vars[sp_name].append(index)
        else:
            sp_vars[sp_name] = []
            sp_vars[sp_name].append(index)
        # genus
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
    fp = open(output_dir + "SH_" + name + ".fas", "w")
    for index in sh_vars[name]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash_md5[index] + "|SampleID_" + samples_mapping[s] + "|sh_" + name + "|marker_" + markers[index] + "|abund_"  + a +"_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()
    zip_file(output_dir + "SH_" + name + ".fas")

print("SH variants were written...")

# sp
for name in sp_vars:
    if "_sp." not in name:
        i = 1
        fp = open(output_dir + "species_" + name + ".fas", "w")
        for index in sp_vars[name]:
            for i in range(len(samples[index])):
                s = samples[index][i]
                a = abundances[index][i]
                fp.write(">" + hash_md5[index] + "|SampleID_" + samples_mapping[s] + "|species_" + name + "|marker_" + markers[index] + "|abund_" + a + "_total_" + str(total[s]) + "\n")
                fp.write(sequences[index] + "\n")
                i += 1
        fp.close()
        zip_file(output_dir + "species_" + name + ".fas")
    else:
        print("Ignored: "+name)

print("Species variants were written...")

# gen
for name in gen_vars:
    i = 1
    fp = open(output_dir + "genus_" + name + ".fas", "w")
    for index in gen_vars[name]:
        for i in range(len(samples[index])):
            s = samples[index][i]
            a = abundances[index][i]
            fp.write(">" + hash_md5[index] + "|SampleID_" + samples_mapping[s] + "|genus_" + name + "|marker_" + markers[index] + "|abund_" + a + "_total_" + str(total[s]) + "\n")
            fp.write(sequences[index] + "\n")
            i += 1
    fp.close()
    zip_file(output_dir + "genus_" + name + ".fas")

print("Genus variants were written...")

########
print("Done :)")


