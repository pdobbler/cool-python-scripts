__author__ = 'Wietrack 2026'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]
blast_out6 = sys.argv[2]
group_taxa = sys.argv[3]
out_table = sys.argv[4]
release = sys.argv[5]               # 2

thr_s = 197.0
thr_g = 190.0

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def convert_cl_code(value):
    # Remove "CL" prefix, then remove leading zeroes
    result = value[2:].lstrip('0')

    # Return "0" if the remaining value was all zeroes
    return result or '0'

#############################################

#groupID taxonomy
#g0000001        d__Bacteria;p__Bacillota_I;c__Bacilli_A;o__Alicyclobacillales;f__Alicyclobacillaceae_368324;g__Alicyclobacillus_A_368256;s__
#g0000002        d__Bacteria;p__Bacillota_A_368345;c__Clostridia_258483;o__Oscillospirales;f__Ruminococcaceae;g__Gemmiger_A_73129;s__

taxonomy = {}
for line in openfile(group_taxa):
    values = line.strip().split("\t")
    taxonomy[values[0]] = values[1]

print("groups taxonomy loaded...")

#############################################

#CL34729200|8fa4d5c68ea4bbcd124d1a64c0907b69;samples=1;relabund_sum=0.0000128559;size=1|100.0    g0083103        91.304  253     22      0       1       253  1 253     2.61e-95        346
#CL34729201|ca0165ab368aad953816a8436f04675c;samples=1;relabund_sum=0.0000128559;size=1|100.0    g0070140        92.857  252     18      0       2       253  2 253     2.00e-101       366

#CL34729200|8fa4d5c68ea4bbcd124d1a64c0907b69  RS-GCF-000455125.1-NZ-AUWR01000057.1|d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacteriumvirginiense 91.304  253     22      0       1       253     517     769     5.76e-94        346
#CL34729201|ca0165ab368aad953816a8436f04675c  MJ034-2-barcode67-umi13745bins-ubs-8|d__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__UBA2199;f__UBA2199;g__UBA2199;s__UBA2199sp002436105     92.857  252     18      0       2       253     520     771     4.42e-100       366

# load blast results...
blast_info = {}
for line in openfile(blast_out6):
    values = line.strip().split("\t")
    # storing - pident | length |  evalue | bitscore
    #blast_info[values[0]] = values[1] + "\t" + values[2] + "\t" + str(int(values[7])-int(values[6])+1) + "\t" + values[10] + "\t" + values[11]
    # new                            group               similarity                 qend            qstart                  evalue              bitscore
    blast_info[values[0]] = taxonomy[values[1]] + "\t" + values[2] + "\t" + str(int(values[7])-int(values[6])+1) + "\t" + values[10] + "\t" + values[11]

# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#############################################

#>CL00000001|d829bee4984f82ffc2453212157caf96;samples=25791;relabund_sum=77.2086052524;size=5579362|100.0
#TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG
#>CL00000002|90aca794c7e30b8a77e87f13ffc9a5cc;samples=21372;relabund_sum=40.4904515450;size=2629276|100.0
#TACGAAGGGG

# output
# 1       GB00000001.2    Bradyrhizobium sp000244915      Bradyrhizobium_503372   100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rhizobiales_505101; f__Xanthobacteraceae; g__Bradyrhizobium_503372; s__Bradyrhizobium sp000244915    6472eb8b1e09f892aca2f23182962903

n = 0
ot = open(out_table, 'w')
titleRead = False
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            slen = len(seq)
            # process title
            tparts = title.split(";")[0].split("|")
            clIndex = convert_cl_code(tparts[0])
            clID = tparts[0].replace("CL", "GB") + "." + release
            md5 = tparts[1]
            # process ident
            gen = '-'
            sp = '-'
            tsim = "NA"
            tcov = "NA"
            taxonomy = "unidentified"
            if blast_info.has_key(title):
                bl_info = blast_info[title].split("\t")
                sim = float(bl_info[1])
                cov = (float(bl_info[2])/slen)*100
                if sim + cov >= thr_g:
                    if len(bl_info[0].split(';')[5])>3:
                        gen = bl_info[0].split(';')[5].replace("g__", "")
                if sim + cov >= thr_s:
                    if len(bl_info[0].split(';')[6])>3:
                        sp = bl_info[0].split(';')[6].replace("s__", "")
                taxonomy = bl_info[0].replace(";", "; ")
                tsim = str(sim)
                tcov = str(cov)
            else:
                n = n + 1
            # write
            ot.write(clIndex+"\t"+clID+"\t"+sp+"\t"+gen+"\t"+tsim+"\t"+tcov+"\t"+taxonomy+"\t"+md5+"\n")
ot.close()

print("Done - taxonomy table created - unidentified: "+str(n))
