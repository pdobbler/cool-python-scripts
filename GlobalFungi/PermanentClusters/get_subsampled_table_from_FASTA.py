__author__ = 'vetrot'

import sys
import random
import gzip

fasta_file = sys.argv[1]
size = int(sys.argv[2])
out_file = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#------------------------------------------

samples = {}
for line in openfile(fasta_file):
    if line[0] == '>':
        sample = line.rstrip().split('|')[1]
        cluster = line.rstrip().split('|')[2]
        if samples.has_key(sample):
            sample_list = samples[sample]
            sample_list.append(cluster)
        else:
            samples[sample] = [cluster]

print("subsample samples to "+str(size))

clusters = {}
sample_names = []
for sample in samples:
    sample_list = samples[sample]
    if len(sample_list) >= size:
        sample_names.append(sample)
        random_list = random.sample(sample_list, size)
        for c in random_list:
            if clusters.has_key(c):
                smap = clusters[c]
                if smap.has_key(sample):
                    smap[sample] += 1
                else:
                    smap[sample] = 1
                clusters[c] = smap
            else:
                smap = {}
                smap[sample] = 1
                clusters[c] = smap

print("samples were subsampled - number of samples: "+str(len(sample_names))+ " clusters "+str(len(clusters)))

fp = open(out_file, 'w')
line = "cluster_name"
for sample in sample_names:
    line += "\t" + sample
fp.write(line + "\n")

for c in clusters:
    line = c
    smap = clusters[c]
    for sample in sample_names:
        if smap.has_key(sample):
            line += "\t" + str(smap[sample])
        else:
            line += "\t0"
    fp.write(line + "\n")
fp.close()

print("DONE :]")

