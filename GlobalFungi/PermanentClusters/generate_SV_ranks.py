__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]

class Variant:
    def __init__(self, sequence, sample, study):
        self.samples  = set()
        self.samples.add(sample)
        self.papers   = set()
        self.papers.add(study)
        self.sequence = sequence
        self.size = 1
        self.rank = 0.0

    def setRank(self, max_var_size, max_sample_size, max_paper_size):
        self.rank = self.getSVRank(max_var_size, max_sample_size, max_paper_size)

    def add(self, sample, study):
        self.samples.add(sample)
        self.papers.add(study)
        self.size += 1

    def getInfo(self, max_var_size, max_sample_size, max_paper_size):
        line = hashlib.md5(self.sequence.encode()).hexdigest()
        line += '|' + "V_" + str(self.size)
        line += '|' + "S_" + str(len(self.samples))
        line += '|' + "P_" + str(len(self.papers))
        line += '|' + "r_" + str(self.rank)
        return line

    def getSeq(self):
        return self.sequence

    def getSize(self):
        return self.size

    def getSampleSize(self):
        return len(self.samples)

    def getPaperSize(self):
        return len(self.papers)

    def getSVRank(self, max_var_size, max_sample_size, max_paper_size):
        # SV Rank = 0.4 * (seq_vars_size/33910054) + 0.5 * (sample_size/22248) + 0.1 * (study_size / 274)
        return 0.4 * (self.size / float(max_var_size)) + 0.5 * (len(self.samples) / float(max_sample_size)) + 0.1 * (len(self.papers) / float(max_paper_size))

    def isQualified(self):
        if len(self.papers)>=5:
            return True
        if len(self.samples)>=10:
            return True
        if self.size>=100 and len(self.samples)>=5:
            return True
        if self.size>=100 and len(self.papers)>=3:
            return True
        return False

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >GF4S05056b|Sun_2021_PK|ERR4885923.1221774
i = 0
sequences = {}
max_sample_size = 0
max_paper_size  = 0
max_var_size = 0
for line in openfile(fasta_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            # process SV
            vals = title.split('|')
            if sequences.has_key(seq):
                v = sequences[seq]
                v.add(vals[0], vals[1])
                #
                if max_var_size < v.getSize():
                    max_var_size = v.getSize()
                if max_sample_size < v.getSampleSize():
                    max_sample_size = v.getSampleSize()
                if max_paper_size < v.getPaperSize():
                    max_paper_size = v.getPaperSize()
                #
                sequences[seq] = v
            else:
                v = Variant(seq, vals[0], vals[1])
                sequences[seq] = v
            i += 1
            #if i>10000:
            #    print("testing version...")
            #    break
print("Sequence processed: "+str(i))

ut = list(sequences.values())
for v in ut:
    v.setRank(max_var_size, max_sample_size, max_paper_size)

print("Ranks were set...")
print("Sorting by rank...")
ut.sort(key=lambda x: x.rank, reverse=True)
print("...done")
print("Writing sorted files...")
sum = 0
q = 0
n = 0
fpQ = open(fasta_file+".qualified", 'w')
fpN = open(fasta_file+".nonqualified", 'w')
for v in ut:
    if v.isQualified():
        fpQ.write(">"+v.getInfo(max_var_size, max_sample_size, max_paper_size)+"\n")
        fpQ.write(v.getSeq() + "\n")
        q += 1
    else:
        fpN.write(">" + v.getInfo(max_var_size, max_sample_size, max_paper_size) + "\n")
        fpN.write(v.getSeq() + "\n")
        n += 1
    sum += v.getSize()
fpN.close()
fpQ.close()

print("Sequence saved: "+str(sum)+" qualified variants: "+str(q)+" non-qualified variant: "+str(n))
print("max_var_size: "+str(max_var_size))
print("max_sample_size: "+str(max_sample_size))
print("max_paper_size: "+str(max_paper_size))

