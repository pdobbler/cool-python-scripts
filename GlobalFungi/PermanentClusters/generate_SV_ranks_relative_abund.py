__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]
split = sys.argv[2].lower() == 'true'  # Convert the string to a boolean

if split == true:
    print("Daset will be split to Qualified and Nonqualified sequences...")
else:
    print("No split applied!")



def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

class Variant:
    def __init__(self, sequence, sample, study):
        self.samples  = {}
        self.samples[sample] = 1
        self.papers   = set()
        self.papers.add(study)
        self.sequence = sequence
        self.size = 1
        self.rank = 0.0

    def add(self, sample, study):
        if self.samples.has_key(sample):
            self.samples[sample] += 1
        else:
            self.samples[sample] = 1
        self.papers.add(study)
        self.size += 1

    def setRank(self, samples_max):
        self.rank = self.getSVRank(samples_max)

    def getInfo(self):
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

    def getSVRank(self, samples_max):
        rel_abund_sum = 0.0
        for sample in self.samples:
            rel_abund_sum += self.samples[sample]/ float(samples_max[sample])
        return rel_abund_sum
        # SV Rank = 0.4 * (seq_vars_size/33910054) + 0.5 * (sample_size/22248) + 0.1 * (study_size / 274)
        #return 0.4 * (self.size / float(max_var_size)) + 0.5 * (len(self.samples) / float(max_sample_size)) + 0.1 * (len(self.papers) / float(max_paper_size))

    def isQualified(self):
        if len(self.samples)>=5:
            return True
        #if len(self.papers)>=5:
        #    return True
        #if self.size>=100 and len(self.samples)>=5:
        #    return True
        #if self.size>=100 and len(self.papers)>=3:
        #    return True
        return False

# >GF4S05056b|Sun_2021_PK|ERR4885923.1221774
i = 0
sequences = {}
samples_max = {}
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

            # count samples sequences
            sampleName = vals[0]
            studyName = vals[1]
            if samples_max.has_key(sampleName):
                samples_max[sampleName] += 1
            else:
                samples_max[sampleName] = 1

            # process variant...
            if sequences.has_key(seq):
                v = sequences[seq]
                v.add(sampleName, studyName)
                sequences[seq] = v
            else:
                sequences[seq] = Variant(seq, sampleName, studyName)
            i += 1
            #if i>10000:
            #    print("testing version...")
            #    break
print("Sequence processed: "+str(i))

ut = list(sequences.values())
for v in ut:
    v.setRank(samples_max)

print("Ranks were set...")
print("Sorting by rank...")
ut.sort(key=lambda x: x.rank, reverse=True)
print("...done")
print("Writing sorted files...")
sum = 0
q = 0
n = 0
if split:
    fpQ = open(fasta_file+".qualified", 'w')
    fpN = open(fasta_file+".nonqualified", 'w')
    for v in ut:
        if v.isQualified():
            fpQ.write(">"+v.getInfo()+"\n")
            fpQ.write(v.getSeq() + "\n")
            q += 1
        else:
            fpN.write(">" + v.getInfo() + "\n")
            fpN.write(v.getSeq() + "\n")
            n += 1
        sum += v.getSize()
    fpN.close()
    fpQ.close()
else:
    fpQ = open(fasta_file+".all", 'w')
    for v in ut:
        fpQ.write(">"+v.getInfo(max_var_size, max_sample_size, max_paper_size)+"\n")
        fpQ.write(v.getSeq() + "\n")
        q += 1
        sum += v.getSize()
    fpQ.close()

print("Sequence saved: "+str(sum)+" qualified variants: "+str(q)+" non-qualified variant: "+str(n))

