__author__ = 'avetrot'

import sys
import gzip

FASTA = sys.argv[1]
motive = sys.argv[2]
mismatches = int(sys.argv[3])

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# loading of barcodes...
print("motive " + motive + " mismatch allowed " + str(mismatches))

def compareStrings(s1, s2):
    mm = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            mm += 1
    return mm

def search(sequence, adapt, mismatch):
    for i in range(len(sequence)):
        part = sequence[i:]
        if len(part) < len(adapt):
            break
        p = part[0:len(adapt)]
        mm = compareStrings(p, adapt)
        if mm <= mismatch and len(p) > (mismatch + 1):
            return i
    return -1

# processing the fasta
filled = False
count = 0
no_hit = 0

# Open output file in gzip mode
with gzip.open(FASTA + ".SECONDARYCUT.fa.gz", 'wt') as fp:
    for n, line in enumerate(openfile(FASTA)):
        if n % 2 == 0:
            r1_0 = line.rstrip()
        else:
            if n % 2 == 1:
                r1_1 = line.rstrip()
                filled = True
        if filled:
            filled = False
            # search for motives
            title = r1_0.split(' ')[0]
            seq = r1_1
            pos = search(seq, motive, mismatches)
            if pos > -1:
                fp.write(title + '|POS=' + str(pos) + '\n')
                fp.write(seq[:pos] + '\n')
                count += 1
            else:
                no_hit += 1

print(FASTA + " Total seqs: " + str(count + no_hit) + " Total found: " + str(count) + " NO HIT: " + str(no_hit))
