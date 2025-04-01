import sys
import random
import gzip
from collections import defaultdict

# Inputs
size = int(sys.argv[1])           # number of sequences per sample
fasta_file = sys.argv[2]          # input fasta file
output_file = sys.argv[3]         # output file for selected sequences
discarded_file = sys.argv[4]      # file to save names of discarded samples

# open gzip files
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)

# Read all sequences and group them by sample
samples = defaultdict(list)
with openfile(fasta_file) as f:
    while True:
        header = f.readline()
        if not header:
            break
        seq = f.readline()
        if not seq:
            break
        if header.startswith('>'):
            sample_name = header.split('|')[0][1:]  # Get 'GB01001571S' from '>GB01001571S|...'
            samples[sample_name].append((header.strip(), seq.strip()))

# Select and save random sequences
discarded_samples = []
with open(output_file, 'w') as out:
    for sample, seqs in samples.items():
        if len(seqs) >= size:
            selected = random.sample(seqs, size)
            for header, sequence in selected:
                out.write(f"{header}\n{sequence}\n")
        else:
            discarded_samples.append(sample)

# Save discarded sample names
with open(discarded_file, 'w') as dfile:
    for sample in discarded_samples:
        dfile.write(sample + "\n")

print(f"Done. Selected sequences written to {output_file}")
print(f"Discarded samples written to {discarded_file}")
