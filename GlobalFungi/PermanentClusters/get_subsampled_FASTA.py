__author__ = 'vetrot'

import sys

fasta_file = sys.argv[1]
size = int(sys.argv[2])
out_file = sys.argv[3]


n = 0
samples = {}
for line in open(fasta_file):
    if line[0] == '>':
        val = line.rstrip().split('|')
        sample = val[1]
        if samples.has_key(sample):
            sample_list = samples[sample]
            sample_list.append(line.rstrip())
        else:
            samples[sample] = [line.rstrip()]

print("subsample samples to "+str(size))

fp = open(out_file, 'w')
for sample in samples:
    random_list = random.sample(l, size)
    for i in random_list:
        fp.write(i + "\n")
fp.close()

print("DONE :]")

