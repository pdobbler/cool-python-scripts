__author__ = 'vetrot'

import sys

original_table = sys.argv[1]        # MT_Bavarian_Ga0566201_TAX_CAZy_KEGG_KOG_RENAMED.tab
tax_column = int(sys.argv[2])       # 18
KOG_column = int(sys.argv[3])       # 24
new_fungal_hits = sys.argv[4]       # ga0566201_proteins_nolinebreaks_JGI_FUN_20240403_best.txt
new_taxonomy = sys.argv[5]          # taxonomy_2024.txt
tax_tree = sys.argv[6]              # TAX_tree.tab

taxa_set = set()
for line in open(original_table):
    vals = line.rstrip().split('\t')
    print("taxonomy column to compare: " + vals[tax_column]+" bitscore> "+vals[tax_column-1])
    print("KOG column to replace: " + vals[KOG_column]+" evalue> "+vals[KOG_column-1])
    break



# jgi|Tubae1|2744|KOG0633|2.6.1.9 41.9
gene_hits = {}
for line in open(new_fungal_hits):
    vals = line.rstrip().split('\t')
    gene_hits[vals[0]] = vals[1].split('|')[0]+'|'+vals[11]+'|'+vals[1].split('|')[1]+'|'+vals[10]
    #print(gene_hits[vals[0]])

print("hits loaded...")

new_tax_tree = {}
for line in open(new_taxonomy):
    vals = line.rstrip().split('\t')
    new_tax_tree[vals[6]] = line.rstrip()

print("new taxonomy loaded...")

taxons_ori = {}
for line in open(tax_tree):
    vals = line.rstrip().split('\t')
    taxons_ori[vals[6]] = 0

print("original tax tree loaded...")

uinfo = open("update_info.txt", "w")
uinfo.write("ori_taxon\tnew_taxon\tori_bitscore\tnew_bitscore\tori_KOG\tnew_KOG\tori_eval\tnew_eval\n")
out_file = open(original_table + ".updated", "w")
for line in open(original_table):
    vals = line.rstrip().split('\t')
    if gene_hits.has_key(vals[0]):
        hit_vals = gene_hits[vals[0]].split('|')
        if float(vals[tax_column-1])<float(hit_vals[1]):
            uinfo.write(vals[tax_column] + "\t" + hit_vals[0]+"\t" + vals[tax_column-1] + "\t" +  hit_vals[1]+"\t"+vals[KOG_column] +"\t["+ hit_vals[2]+"]\t"+vals[KOG_column-1] + "\t" + hit_vals[3]+ "\n")
            #print("replacing original value: " + vals[tax_column] + " by " + hit_vals[0]+" bitscore " + vals[tax_column-1] + " vs " +  hit_vals[1]+" | KOG> "+vals[KOG_column] +" vs ["+ hit_vals[2]+"] evalue "+vals[KOG_column-1] + " vs " + hit_vals[3])
            vals[tax_column] = hit_vals[0]
            vals[tax_column-1] = hit_vals[1]
            vals[KOG_column] = "[" + hit_vals[2] + "]"
            vals[KOG_column-1] = hit_vals[3]
            # update tax tree...
            if taxons_ori.has_key(hit_vals[0]):
                taxons_ori[hit_vals[0]] = 2
            else:
                taxons_ori[hit_vals[0]] = 1
        else:
            taxons_ori[vals[tax_column]] = 2
    # write to a new table...
    out_file.write('\t'.join(vals) + "\n")
out_file.close()
uinfo.close()

print("updated table was saved...")

nt = open("TAX_tree_updated.tab", "w")
for n, line in enumerate(open(tax_tree)):
    if n == 0:
        nt.write(line.rstrip()+"\n")
    else:
        vals = line.rstrip().split('\t')
        if taxons_ori[vals[6]] > 0:
            nt.write(line.rstrip()+"\n")
# add new taxons...
for taxon in taxons_ori:
    if taxons_ori[taxon] == 1:
        nt.write(new_tax_tree[taxon]+"\n")
nt.close()

print("DONE :)")




