__author__ = 'vetrot'

import sys
import os
import gzip

taxonomy_ident = sys.argv[1]        # CLUSTERS_INFO_TOTAL_FINAL.txt

thr_s = 197.0
thr_g = 190.0

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# clName  md5     totSeqSize      totRepSize      HIT     TAXONOMY        SIMILARITY      COVERAGE        EVALUE  BITSCORE        SEQ     totalSeqs       totalSamples    totalStudies    repSamples      repStudies      clusterRelAbund_sum     clusterVariants clusterVariantSize_sum
# CL00000001      6472eb8b1e09f892aca2f23182962903        27490562        6712798 GB-GCA-004799405.1-SSMW01000135.1       d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales_505101;f__Xanthobacteraceae;g__Bradyrhizobium_503372;s__Bradyrhizobiumsp000244915   100.000 100.0   1.18e-130       468     TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG   27490562        28379   169     23821   164     448.303121968   2994881 27490562
# CL00000002      90aca794c7e30b8a77e87f13ffc9a5cc        27542758        2629276 MJ034-1-barcode32-umi22187bins-ubs-8    d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales_505101;f__Xanthobacteraceae;g__Variibacter;s__      100.000 100.0   1.18e-130       468     TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGCTTTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGAAGCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG   27542758        27557   166     21372   155     421.793604574   2354933 27542758

sp_taxonomy = {}
gen_taxonomy = {}
i = 0
fp = open("TAXONOMY_CLUSTERS.txt", 'w')
for line in openfile(taxonomy_ident):
    line = line.rstrip()
    if i > 0:
        parts = line.split("\t")
        taxons = parts[5]
        md5 = parts[1]
        cl = parts[0]
        genus = '-'
        species = '-'
        if taxons == 'UNKNOWN':
            taxons = 'unidentified'
        else:
            sim = float(parts[6])
            cov = float(parts[7])
            if sim + cov >= thr_g:
                gen = taxons.split(';')[5].replace("g__", "")
                if not gen == '':
                    genus = gen
                    gen_taxonomy[genus] = taxons
                if sim + cov >= thr_s:
                    sp = taxons.split(';')[6].replace("s__", "")
                    if not sp == '':
                        species = sp
                        sp_taxonomy[species] = taxons
        fp.write(cl + '\t' + species + '\t' + genus + '\t' + str(sim) + '\t' + str(cov) + '\t' + taxons.replace(";", "; ") + '\t' + md5 + '\n')
    i += 1
fp.close()

print("Taxonomy loaded: " + str(i) + " - the size should be the same as clusters!")

fp = open("TAXONOMY_SP_GEN.txt", 'w')
for species in sp_taxonomy:
    taxons = sp_taxonomy[species]
    gen = taxons.split(';')[5].replace("g__", "")
    if gen_taxonomy.has_key(gen):
        del gen_taxonomy[gen]
    fp.write('sp\t' + species + '\t' + taxons + '\n')

for genus in gen_taxonomy:
    taxons = gen_taxonomy[genus]
    fp.write('gen\t' + genus + '\t' + taxons + '\n')
    
fp.close()

print("Done :]")
