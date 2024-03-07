__author__ = 'vetrot'

import sys
import datetime
from Bio import Align

fasta_file = sys.argv[1]
sim_treshold = float(sys.argv[2])

print("Grouping based on similarity threshold: " + str(sim_treshold))
t1 = datetime.datetime.now()

aligner = Align.PairwiseAligner()
aligner.mode = "global"

print("Pairwise Aligner mode: " + aligner.mode)

class Sequence:
    def __init__(self, title, sequence, cluster, seed):
        self.title = title
        self.sequence = sequence
        self.cluster  = cluster
        self.seed = seed
        self.sim = 0.0
        self.length = len(sequence)

    def getCluster(self):
        return self.cluster

    def getTitle(self):
        return self.title

    def getSeq(self):
        return self.sequence

    def getLen(self):
        return self.length

    def getSim(self):
        return self.sim

    def isSeed(self):
        return self.seed

    def setCluster(self, cluster, seed, sim):
        self.cluster = cluster
        self.seed = seed
        self.sim = sim


seeds = {}
sequences = []
for line in open(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            sequences.append(Sequence(title, seq, "", False))

t2 = datetime.datetime.now()
dif = t2-t1
print("All sequences loaded "+str(len(sequences))+" total time "+str(dif.total_seconds())+" sec")
t1 = datetime.datetime.now()

#
seeds = []
cluster = 1
s = sequences[0]
s.setCluster(cluster,True, 100.0)
seeds.append(s) # first seed
sequences[0] = s
print("Cluster created " + str(cluster) + " - seqs remains: " + str(len(sequences)-1))
#
for i in range(1, len(sequences)):
    s = sequences[i]
    found = False
    for seed in seeds:
        dl = float(abs(seed.getLen() - s.getLen()))
        lensim = 100.0 - (dl / seed.getLen() * 100.0)
        if lensim >= sim_treshold:
            alignments = aligner.align(seed.getSeq(), s.getSeq())
            score = alignments[0].score
            similarity = score / float(seed.getLen()) * 100.0
            if similarity >= sim_treshold:
                s.setCluster(seed.getCluster(), False, similarity)
                found = True
                break
    if not found:
        cluster += 1
        s.setCluster(cluster, True, 100.0)
        seeds.append(s)
        print("Cluster created " + str(cluster) + " - seqs remains: " + str(len(sequences)-(i+1)))
#
t2 = datetime.datetime.now()
dif = t2-t1
print("Done - # clusters "+str(cluster)+" total time "+str(dif.total_seconds())+" sec")

fpC = open(fasta_file + ".clustered2", 'w')
fpS = open(fasta_file + ".seeds2", 'w')
name_max = len(str(cluster))
for s in sequences:
    cname = str(s.getCluster())
    for i in range(name_max-len(cname)):
        cname = "0"+ cname
    cname = "CL" + cname
    title = cname + '|' + s.getTitle()
    seq = s.getSeq()
    #print(cname)
    if s.isSeed():
        fpS.write('>'+title+'|SEED\n')
        fpS.write(seq+'\n')
    fpC.write('>'+title+'|'+str(s.getSim())+'\n')
    fpC.write(seq+'\n')
fpC.close()
fpS.close()

print("Saved :)")

