import sys

d_gene = {}

for f in sys.argv[1:]:
    for line in open(f):
        ch = line[0]
        if ch != '@':
            mg_id = f.split('.txt')[0]
            gene_name = line.rstrip().split('\t')[0]
            gene_length = line.rstrip().split('\t')[1]
            gene = gene_name+"\t"+gene_length
            count = int(line.rstrip().split('\t')[2]) #mapped
            #count = int(dat[3]) #unmapped

            if d_gene.has_key(gene):
                d_gene[gene][mg_id] = count
            else:
                d_gene[gene] = {}
                d_gene[gene][mg_id] = count

fp = open('summary-count-mapped.tsv', 'w')

sorted_samples = sys.argv[1:]

fp.write('ctg_name\tctg_length')

for x in sorted_samples:
    fp.write('\t%s' % x.split('.')[0])

fp.write('\n')

for gene in d_gene:
    fp.write('%s\t' % gene)
    for x in sorted_samples:
        x1 = x.split('.txt')[0]
        if d_gene[gene].has_key(x1):
            fp.write('%s\t' % d_gene[gene][x1])
        else:
            fp.write('0\t')
    fp.write('\n')




