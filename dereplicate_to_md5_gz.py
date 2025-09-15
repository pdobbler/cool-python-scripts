__author__ = 'Wietrack 2016'

import sys
import operator
import hashlib
import gzip

fasta_file = sys.argv[1]
out_fasta = sys.argv[2]
out_table = sys.argv[3]


def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

fasta = {}
loaded = False
titleRead = False

i = 0
n = 0
#load fasta sequences
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            if len(seq)>0:
                fasta[title] = seq
                i=i+1
            else:
                n = n + 1

print str(i)+" sequences loaded correctly - "+str(n)+" sequnces are empty - Omitted!"

#sort the dictionary
sorted_fasta = sorted(fasta.items(), key=operator.itemgetter(1))
print str(len(sorted_fasta))+" sequence variants sorted ..."

#dereplication...
of = open(out_fasta, 'w')
ot = open(out_table, 'w')

n = 1
g = 0
gname=''
derSeq=''
seq_last = ''
seq = ''
for f in sorted_fasta:
    seq = f[1]
    if seq == seq_last:
        n=n+1
    else:
        if derSeq<>'':
            of.write('>'+gname+'\n')
            of.write(derSeq+'\n')
            g = g + 1
        gname = hashlib.md5(seq.encode()).hexdigest()
        derSeq=seq
        n=1
    ot.write(gname+'\t'+f[0]+'\n')
    seq_last = seq
#last sequence...
if derSeq<>'':
            of.write('>'+gname+'\n')
            of.write(derSeq+'\n')
            g = g + 1

of.close()
ot.close()

print 'Dereplication is done - '+str(g)+' groups from '+str(len(sorted_fasta))+' seqs'



