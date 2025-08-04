__author__ = 'vetrot'

import sys
import os
import gzip

cluster_ident = sys.argv[1]         # VARIANTS_TABLE_CLUSTER_PAIRS.txt
taxonomy_ident = sys.argv[2]        # SEEDS_97.0_WORKING_NAMES_TAXONOMY.txt


thr_s = 197.0
thr_g = 190.0


def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# GB00000040.1    p__Actinomycetota|CL000009|087e50a409374f400614f96d9c161f62|V_513560|S_230|P_6|r_7.57535829642|SEED
# 1       GB00000001.1    6472eb8b1e09f892aca2f23182962903
md5_clusters = {}
i = 0
for line in openfile(cluster_ident):
    parts = line.rstrip().split('\t')
    md5_clusters[parts[2]] = parts[0] + '\t' + parts[1]
    i += 1

print("Clusters loaded: " + str(i))

# SH0920751.10FU  Fungi   Basidiomycota   Agaricomycetes  Thelephorales   Thelephoraceae  Tomentella      Tomentella sp.  12147
# 020a94d69ebfddc17d8e0b45d9d140ee        CL000027_p__Actinomycetota      1024746 3497    67      3.29711953244   100.000 100.0   d__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Streptomycetales_400645; f__Streptomycetaceae_400641; g__Streptomyces_400150; s__Streptomyces africanus
# dc6a82a18445efcd5facad6ab0c8812a        CL00008_p__Eremiobacterota      366771  222     6       0.659359130686  99.209  100.0   d__Bacteria; p__Eremiobacterota; c__Xenobia; o__Xenobiales; f__; g__; s__
# 6eb9daca50e4b0a2ddc361392e1f93af        CL20407_NO_HIT  2       1       1       1.53563832646e-05       -       -       NO_HIT
sp_taxonomy = {}
gen_taxonomy = {}
i = 0
fp = open("TAXONOMY_CLUSTERS.txt", 'w')
for line in openfile(taxonomy_ident):
    line = line.rstrip()
    parts = line.split("\t")
    taxons = parts[8]
    md5 = parts[0]
    if md5_clusters.has_key(md5):
        cl = md5_clusters[md5]
        genus = '-'
        species = '-'
        if taxons == 'NO_HIT':
            taxons = 'unidentified'
        else:
            sim = float(parts[6])
            cov = float(parts[7])
            if sim + cov >= thr_g:
                gen = taxons.split(';')[5].replace(" g__", "")
                if not gen == '':
                    genus = gen
                    gen_taxonomy[genus] = taxons
                if sim + cov >= thr_g:
                    sp = taxons.split(';')[6].replace(" s__", "")
                    if not sp == '':
                        species = sp
                        sp_taxonomy[species] = taxons
        fp.write(cl + '\t' + species + '\t' + genus + '\t' + str(sim) + '\t' + str(cov) + '\t' + taxons + '\t' + md5 + '\n')
        i += 1
fp.close()

print("Taxonomy loaded: " + str(i) + " - the size should be the same as clusters!")

fp = open("TAXONOMY_SP_GEN.txt", 'w')
for species in sp_taxonomy:
    taxons = sp_taxonomy[species]
    gen = taxons.split(';')[5].replace(" g__", "")
    if gen_taxonomy.has_key(gen):
        del gen_taxonomy[gen]
    fp.write('sp\t' + species + '\t' + taxons + '\n')

for genus in gen_taxonomy:
    taxons = gen_taxonomy[genus]
    fp.write('gen\t' + genus + '\t' + taxons + '\n')
    
fp.close()

print("Done :]")
