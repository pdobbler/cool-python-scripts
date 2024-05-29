__author__ = 'vetrot'

import sys
import datetime
from Bio import Align
from itertools import product  # Importing product from itertools

fasta_file = sys.argv[1]
sim_treshold = float(sys.argv[2])

score_thr = (100.0-sim_treshold)/100.0

# Convert command line string argument to boolean
if sys.argv[3].lower() in ['true', '1']:
    show_alignments = True
elif sys.argv[3].lower() in ['false', '0']:
    show_alignments = False
else:
    raise ValueError("Invalid boolean value. Use 'true'/'1' or 'false'/'0'.")

print("Grouping based on similarity threshold: " + str(sim_treshold))
print("Show alignments: " + str(show_alignments))
t1 = datetime.datetime.now()

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
#aligner.match_score = 1
#aligner.mismatch_score = -1


print("Pairwise Aligner mode: " + aligner.mode)
print("Pairwise Aligner open_gap_score: " + str(aligner.open_gap_score))
print("Pairwise Aligner extend_gap_score: " + str(aligner.extend_gap_score))
print("Pairwise Aligner match_score: " + str(aligner.match_score))
print("Pairwise Aligner mismatch_score: " + str(aligner.mismatch_score))

def generate_all_combinations(characters, length):
    return [''.join(p) for p in product(characters, repeat=length)]

def get_kmer_default_lib():
    combinations = generate_all_combinations('ACTG', 4)
    all_kmers = {}
    # Výpis všech kombinací
    for combo in combinations:
        all_kmers[combo] = 0
    return all_kmers  # Return the dictionary

def generate_kmers(string):
    kmers = get_kmer_default_lib()
    for i in range(len(string) - 3):
        kmer = string[i:i+4]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers  # Return the kmer dictionary

def compare_profiles(p1, p2, words):
    seq_score = 0
    for word in words:
        w1 = p1.get(word, 0)  # Use get method to handle missing keys
        w2 = p2.get(word, 0)  # Use get method to handle missing keys
        seq_score += abs(w2 - w1)
    return seq_score

class Sequence:
    def __init__(self, title, sequence, cluster, seed):
        self.title = title
        self.sequence = sequence
        self.cluster  = cluster
        self.seed = seed
        self.sim = 0.0
        self.length = len(sequence)
        self.kmers = generate_kmers(sequence)

    def getProfile(self):
        return self.kmers

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

words = generate_all_combinations('ACTG', 4)
#
seeds = []
cluster = 1
s = sequences[0]
s.setCluster(cluster,True, 100.0)
seeds.append(s) # first seed
sequences[0] = s
print("Cluster created " + str(cluster) + " - seqs remains: " + str(len(sequences)-1) + " seed name " + s.title)
#
for i in range(1, len(sequences)):
    s = sequences[i]
    found = False
    if show_alignments:
        print("")
        print("Testing "+s.title)
        print("")
    for seed in seeds:
        dl = float(abs(seed.getLen() - s.getLen()))
        lensim = 100.0 - (dl / seed.getLen() * 100.0)
        if lensim >= sim_treshold:
            # check kmers
            score = (compare_profiles(seed.getProfile(), s.getProfile(), words)/9)/seed.getLen()
            #print("score: "+str(score) + " score thr: " +str(score_thr))
            if score < score_thr:
                alignments = aligner.align(seed.getSeq(), s.getSeq())
                alignment_str = str(alignments[0])
                match_line = alignment_str.split('\n')[1]
                mm = len(match_line)-match_line.count('|')
                similarity = float(len(match_line)-mm) / float(len(match_line)) * 100.0
                # show alignments...
                if show_alignments:
                    print("SEED: "+seed.title+" vs. "+s.title+" seed-len: "+str(seed.getLen())+" s-len: "+str(s.getLen())+" score "+str(alignments[0].score)+" sim: "+str(similarity)+" align.len: "+str(len(match_line))+" identical: "+str(match_line.count('|'))+" mismatch: "+str(len(match_line)-match_line.count('|')))
                    print(alignments[0])
                # cluster found...
                if similarity >= sim_treshold:
                    s.setCluster(seed.getCluster(), False, similarity)
                    found = True
                    break
        else:
            if show_alignments:
                print("SEED: "+seed.title+" vs. "+s.title+" seed-len: "+str(seed.getLen())+" s-len: "+str(s.getLen()) + " >>> omitting alignment - similarity based on length: " + str(lensim))
                print("")
    if not found:
        cluster += 1
        s.setCluster(cluster, True, 100.0)
        seeds.append(s)
        print("Cluster created " + str(cluster) + " - seqs remains: " + str(len(sequences)-(i+1))+" seed name "+s.title)

#
t2 = datetime.datetime.now()
dif = t2-t1
print("Done - # clusters "+str(cluster)+" total time "+str(dif.total_seconds())+" sec")

fpC = open(fasta_file + ".clustered3", 'w')
fpS = open(fasta_file + ".seeds3", 'w')
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

