__author__ = 'vetrot'

import sys

in_file = sys.argv[1]
read_size = int(sys.argv[2])
out_file = sys.argv[3]

fp = open(out_file, 'w')
for n, line in enumerate(open(in_file)):
    if n>0:
        gene_name = line.rstrip().split('\t')[0]
        gene_length = int(line.rstrip().split('\t')[1])
        if gene_length>0:
            new_line = gene_name+'\t'+line.rstrip().split('\t')[1]
            for x in range(2, len(line.rstrip().split('\t'))):
                reads_count = float(line.rstrip().split('\t')[x])
                norm_val = (read_size * reads_count)/gene_length
                #print gene_name+" "+str(x)+"   "+line.rstrip().split('\t')[x]+"   %.5f" %(norm_val)
                new_line = new_line +'\t'+ str(norm_val)
            fp.write('%s\n' % new_line)
        else:
            print "WARNING: gene length is 0 bp - "+gene_name
    else:
        fp.write('%s\n' % line.rstrip())
print "done..."
fp.close()



