__author__ = 'vetrot'

import sys
import operator
import collections
import gzip

fasta_file = sys.argv[1]
no_singletons = sys.argv[2]
selected_otus = sys.argv[3] # or '-' to ignore


def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


nosingle = False
if no_singletons.lower() == 'true':
    nosingle = True
if nosingle == True:
    print('output will ignore singletons')
else:
    print('output will not ignore singletons')

otus_list = {}
if selected_otus == "-":
    print('all OTUs will be used')
else:
    for line in open(selected_otus):
        otus_list[line.rstrip()] = 0
    print("OTU list size: "+str(len(otus_list)))

# >GF4S04647b|Sun_2021_PK|ERR4885514.335925.|.|.|OTU0000002

# fill table...
tab_dict_otus = {}
tab_dict_samples = {}
for line in openfile(fasta_file):
    if line[0] == ">":
        vals = line[1:].rstrip().split('|')
        # sum otus...
        otu_name = vals[len(vals)-1]
        if otus_list.has_key(otu_name) or not otus_list:
            # count otu hits...
            if not otus_list:
                if otus_list.has_key(otu_name):
                    otus_list[otu_name] += 1
                else:
                    otus_list[otu_name] = 0
            else:
                otus_list[otu_name] += 1
            #print("OTU: "+otu_name)
            # sum samples...
            sample_name = vals[0]
            #print("Name: "+sample_name)
            # fill otu dict...
            if tab_dict_otus.has_key(otu_name):
                samples_counts = tab_dict_otus[otu_name]
                if samples_counts.has_key(sample_name):
                    samples_counts[sample_name] = samples_counts[sample_name] + 1
                else:
                    samples_counts[sample_name] = 1
                tab_dict_otus[otu_name] = samples_counts
            else:
                samples_counts = {sample_name: 1}
                tab_dict_otus[otu_name] = samples_counts
            # fill otu dict...
            if tab_dict_samples.has_key(sample_name):
                otus_counts = tab_dict_samples[sample_name]
                if  otus_counts.has_key(otu_name):
                    otus_counts[otu_name] =  otus_counts[otu_name] + 1
                else:
                    otus_counts[otu_name] = 1
                tab_dict_samples[sample_name] = otus_counts
            else:
                otus_counts = {otu_name: 1}
                tab_dict_samples[sample_name] = otus_counts

# write table...
fp = open("commpressed_otu_tab.otus_to_samples.nosingle_"+str(nosingle)+".txt", "w")
fp.write("OTU\tsample\tabundance\n")
sum_tot = 0
single_otus = {}
for otu_name in tab_dict_otus:
    samples_counts = tab_dict_otus[otu_name]
    ss = []
    aa = []
    sum = 0
    for sample_name in samples_counts:
        ss.append(sample_name)
        aa.append(str(samples_counts[sample_name]))
        sum += samples_counts[sample_name]
    sum_tot += sum
    if nosingle:
        if sum>1:
            fp.write(otu_name + "\t" + ";".join(ss) + "\t" + ";".join(aa) + "\n")
        else:
            #print(otu_name+" is singleton "+";".join(ss)+" "+";".join(aa))
            single_otus[otu_name] = 0
    else:
        fp.write(otu_name + "\t" + ";".join(ss) + "\t" + ";".join(aa) + "\n")
fp.close()

print("...processed seqs "+str(sum_tot)+" global singletons: "+str(len(single_otus)))

# write table...
fp = open("commpressed_otu_tab.samples_to_otus.nosingle_"+str(nosingle)+".txt", "w")
fp.write("OTU\tsample\tabundance\n")
for sample_name in tab_dict_samples:
    otus_counts = tab_dict_samples[sample_name]
    oo = []
    bb = []
    for otu_name in otus_counts:
        if not single_otus.has_key(otu_name):
            oo.append(otu_name)
            bb.append(str(otus_counts[otu_name]))
    if len(oo)>0:
        fp.write(sample_name + "\t" + ";".join(oo) + "\t" + ";".join(bb) + "\n")
fp.close()

print("Done :)")

for otu_name in otus_list:
    if otus_list[otu_name] == 0:
        print("Warning desired otu "+otu_name+" was not found!!!")

print("Done :)")







