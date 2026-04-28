__author__ = 'vetrot'

import sys
import operator
import collections
import gzip

cl_vars_file = sys.argv[1] # GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz
fasta_file = sys.argv[2]   # GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

##########################
# READ VARIANTS CLUSTERS #
##########################
# >CL00000001|d829bee4984f82ffc2453212157caf96;samples=25791;relabund_sum=77.2086052524;size=5579362|100.0
# TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG
# >CL00000001|6472eb8b1e09f892aca2f23182962903;samples=23821;relabund_sum=111.2113833946;size=6712798|99.60474308300395
# TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG
# >CL00000002|90aca794c7e30b8a77e87f13ffc9a5cc;samples=21372;relabund_sum=40.4904515450;size=2629276|100.0
# TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGCTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGAAGCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG

cl_var_seqs = {}
titleRead = False
for line in openfile(cl_vars_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            title_parts = title.split('|')
            clName = title_parts[0]
            seq = line.strip()
            cl_var_seqs[seq] = clName

print('Clustered variants loaded...')

##########################
# READ SEQS WITH SAMPLES #
##########################
# >GB01020442S|An_2019_1acp_Bact|SRR5920425.1|POS=5|POS=253
# TACAGAGGTCCCGAGCGTTGTTCGGATTCACTGGGCGTAAAGGGTGCGTAGGTGGTGGGGTAAGTCGGATGTGAAATCTCGGAGCTCAACTCCGAAATGGCATTGGAAACTGCCCTGCTAGAGGGTCGGAGGGGGGACTGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACCAGTGGCGAAGGCGAGTCCCTGGACGACACCTGACACTGAGGCACGAAAGCTAGGGGAGCAAACAGG
# >GB01020442S|An_2019_1acp_Bact|SRR5920425.2|POS=7|POS=253
# TACGTAGGGGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAATCTTGCGGCTTAACCGCAAGCGGTCATTGGAAACTGGAGGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGTAGCGAACAGG
# >GB01020442S|An_2019_1acp_Bact|SRR5920425.3|POS=5|POS=253
# AACGTAGGAGGCGAGCGTTATCCGGATTTACTGGGCGTAAAGCGCGTGTAGGCGGTTCGGCAAGTTGGATGTGAAATCTCCTGGCTCAACTGGGAGGGGTCGTTCAATACTACCGGACTTGAGGACAGTAGAGGAAGGTGGAATTCCCGGTGTAGTGGTGAAATGCGTAGATATCGGGAGGAACACCAGTGGCGAAAGCGGCCTTCTGGACTGTTCCTGACGCTCAGACGCGAAAGCTAGGGTAGCAAACGGG

tab_dict_otus = {}
tab_dict_samples = {}
titleRead = False
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            title_parts = title.split('|')
            sample_name = title_parts[0]
            seq = line.strip()
            otu_name = cl_var_seqs[seq]
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

print('Dictionaries filled...')

################
# WRITE TABLES #
################

# write table...
fp = open("compressed_otu_tab_otus_to_samples.txt", "w")
fp.write("OTU\tsample\tabundance\n")
sum_tot = 0
single_otus = {}
for otu_name in tab_dict_otus:
    samples_counts = tab_dict_otus[otu_name]
    ss = []
    aa = []
    otu_sum = 0
    for sample_name in samples_counts:
        ss.append(sample_name)
        aa.append(str(samples_counts[sample_name]))
        otu_sum += samples_counts[sample_name]
    sum_tot += otu_sum
    fp.write(otu_name + "\t" + ";".join(ss) + "\t" + ";".join(aa) + "\n")
fp.close()

print("Done :)")






