__author__ = 'Wietrack 2019'

import sys
import operator
import hashlib
import gzip

fasta_file = sys.argv[1]

################################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
################################################

i = 0
titleRead = False
seq_dict = {}
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
    else:
        if titleRead:
            i = i + 1
            titleRead = False
            seq = line.strip()
            if seq_dict.has_key(seq):
                seq_dict[seq] = seq_dict[seq] + 1
            else:
                seq_dict[seq] = 1

print("Processed sequences " + str(i)+ " to variants " + str(len(seq_dict)))

# Sort by size (value), descending
sorted_seqs = sorted(seq_dict.items(), key=lambda item: item[1], reverse=True)

print("Variants sorted...")

md5_titles = {}
out_file = gzip.open(fasta_file + "_variants.fa.gz", 'wt')
for seq, size in sorted_seqs:
    md5_title = hashlib.md5(seq.encode()).hexdigest()
    out_file.write(">" + md5_title + ";size=" + str(size) + '\n')
    out_file.write(seq +'\n')
    if md5_titles.has_key(md5_title):
        print("Error: nonunique md5 " + md5_title + " seq: " + seq + " size " + str(size)+" !!!")
        md5_titles[md5_title] = md5_titles[md5_title] + 1
    else:
        md5_titles[md5_title] = 1

out_file.close()

print("Done :)")
