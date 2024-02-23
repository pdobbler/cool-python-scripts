import sys

koala_ko = sys.argv[1]  # raw koala output
ko_list = sys.argv[2]   # tresholds file
kegg_multi_out = sys.argv[3]

tresholds = {}
for line in open(ko_list):
    dat = line.rstrip().split('\t')
    tresholds[dat[0]] = dat[1]+"|"+dat[2]

print("thresholds loaded")

gene_dict = {}
for line in open(koala_ko):
    if not line[0] == "#":
        dat = line.rstrip().split()
        gene = dat[0]
        ko_name = dat[2]
        if gene_dict.has_key(gene):
            ko_dict = gene_dict[gene]
        else:
            ko_dict = {}
        #
        if tresholds.has_key(ko_name):
            vars = tresholds[ko_name].split('|')
            if vars[1] == "domain":
                    th = float(dat[8])
                    if th >= float(vars[0]) and float(dat[7]) <= 1.0e-5:
                        # ok save it...
                        eval = dat[7]
                        ko_dict[ko_name] = float(eval)
            if vars[1] == 'full':
                    th = float(dat[5])
                    if th >= float(vars[0]) and float(dat[4]) <= 1.0e-5:
                        # ok save it...
                        eval = dat[4]
                        ko_dict[ko_name] = float(eval)
            gene_dict[gene] = ko_dict
        else:
            print("Fatal error: ko name not in the shitty table")

print("raw table loaded")

fp = open(kegg_multi_out, 'w')
for gene in gene_dict:
    ko_dict = gene_dict[gene]
    if len(ko_dict) > 0:
        sort_orders = sorted(ko_dict.items(), key=lambda x: x[1])
        kos_list = []
        for i in sort_orders:
            #print(i[0], i[1])
            kos_list.append('['+i[0]+']')
        #print("RESULT: "+str(sort_orders[0][1]))
        #print(';'.join(kos_list))
        fp.write(gene + '\t' + str(sort_orders[0][1]) + '\t' + ';'.join(kos_list) + '\n')
fp.close()
print("Done")

