__author__ = 'avetrot'

import sys
import gzip

FASTA_FILE = sys.argv[1]
COMPLETE_BLAST = sys.argv[2]
threshold = float(sys.argv[3]) #sim + cov

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

lengths = {}
titleRead = False
for line in openfile(FASTA_FILE, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            cl_name = title.split('|')[0]
            lengths[cl_name] = len(seq)

print("Sequence lengths loaded... "+str(len(lengths)))


clusters = {}
for line in openfile(COMPLETE_BLAST):
    values = line.rstrip().split('\t')
    #make it smaller
    cl_name = values[0].split('|')[0]
    sh_name = values[1].split('|')[0]
    #
    sim = float(values[2])
    cov = (float(int(values[7]) - int(values[6]) + 1)/float(lengths[cl_name]))*100.0
    if (sim+cov) >= threshold:
        if clusters.has_key(cl_name):
            shs = clusters[cl_name]
            shs.add(sh_name)
            clusters[cl_name] = shs
        else:
            shs = set()
            shs.add(sh_name)
            clusters[cl_name] = shs

print("Blast processed...")

fp = open(FASTA_FILE+".shs", 'w')
fp.write('CLUSTER\tSHs\tsize\tseq_len\n')
for cl_name in clusters:
    shs = clusters[cl_name]
    fp.write(cl_name + '\t' + ';'.join(shs) + '\t' + str(len(shs)) +'\t' +str(lengths[cl_name])+ '\n')
fp.close()

print('DONE :]')

