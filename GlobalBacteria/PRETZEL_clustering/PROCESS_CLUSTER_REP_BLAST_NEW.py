__author__ = 'Wietrack 2019'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]
blast_out6 = sys.argv[2]
out_table = sys.argv[3]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


#############################################

#CL34729200|8fa4d5c68ea4bbcd124d1a64c0907b69  RS-GCF-000455125.1-NZ-AUWR01000057.1|d__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacteriumvirginiense 91.304  253     22      0       1       253     517     769     5.76e-94        346
#CL34729201|ca0165ab368aad953816a8436f04675c  MJ034-2-barcode67-umi13745bins-ubs-8|d__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__UBA2199;f__UBA2199;g__UBA2199;s__UBA2199sp002436105     92.857  252     18      0       2       253     520     771     4.42e-100       366

# load blast results...
blast_info = {}
for line in openfile(blast_out6):
    values = line.strip().split("\t")
    # storing - pident | length |  evalue | bitscore
    blast_info[values[0]] = values[1] + "\t" + values[2] + "\t" + str(int(values[7])-int(values[6])+1) + "\t" + values[10] + "\t" + values[11]

#############################################

#>CL00000001|d829bee4984f82ffc2453212157caf96;samples=25791;relabund_sum=77.2086052524;size=5579362|100.0
#TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG
#>CL00000002|90aca794c7e30b8a77e87f13ffc9a5cc;samples=21372;relabund_sum=40.4904515450;size=2629276|100.0
#TACGAAGGGG

# loop over fasta...
ot = open(out_table, 'w')
loaded = False
titleRead = False
i = 0
n = 0
f = 0
#load fasta sequences
ot.write("clName\tmd5\trepSamples\trepRelAbund_sum\trepSeqs\tHIT\tTAXONOMY\tSIMILARITY\tCOVERAGE\tEVALUE\tBITSCORE\tSEQ\n")
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip().split(";")[0]
        ##########################################
        title_parts = line[1:].strip().split("|")
        tCl = title_parts[0]
        info_parts = title_parts[1].split(";")
        tHex = info_parts[0]
        tRepSample = info_parts[1].split("=")[1]
        tRepRelatA = info_parts[2].split("=")[1]
        tRepSize = info_parts[3].split("=")[1]
        clInfo = tCl + "\t" + tHex + "\t" + tRepSample + "\t" + tRepRelatA + "\t" + tRepSize
        ##########################################
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            slen = len(seq)
            if slen>0:
                if blast_info.has_key(title):
                    values = blast_info[title].split("\t")
                    idents = values[0].split("|")
                    if len(idents)>1:
                        ot.write(clInfo+"\t"+idents[0]+"\t"+idents[1]+"\t"+values[1]+"\t"+str((float(values[2])/slen)*100)+"\t"+ values[3]+"\t"+ values[4]+"\t"+seq+"\n")
                    else:
                        print('Error:' + blast_info[title])
                    f = f + 1
                else:
                    ot.write(clInfo+"\tNO_HIT\tUNKNOWN\t-\t-\t-\t-\t"+seq+"\n")
                i=i+1
            else:
                n = n + 1
ot.close()


print(str(i)+" sequences loaded correctly - "+str(n)+" sequnces are empty - Omitted!")
print("FOUND BLAST RESULTS FOR "+str(f)+" SEQUENCES OUT OF "+str(i))

