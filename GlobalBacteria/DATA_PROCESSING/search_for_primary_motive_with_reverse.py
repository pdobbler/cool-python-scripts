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

def compl_nucl(x):
    return {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }.get(x, 'N')    # 9 is default if x not found

def make_revcompl(seq):
    rev_compl_seq = ''
    for i in range(0,len(seq)):
        rev_compl_seq=rev_compl_seq+compl_nucl(seq[len(seq)-(i+1)].upper())
    return rev_compl_seq

# loading of barcodes...
print("motive " + motive + " mismatch allowed " + str(mismatches))

def compareStrings(s1, s2):
    mm = 0
    for i in range(0, len(s1)):
        if not s1[i] == s2[i]:
            mm += 1
    return mm

def search(sequence, adapt, mismatch):
    for i in range(0, len(sequence)):
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
fwd_count = 0
rev_count = 0
no_hit = 0
fp = open(FASTA + ".PRIMARY.fa", 'w')
for n, line in enumerate(open(FASTA)):
    if n % 2 == 0:
        r1_0 = line.rstrip()
    else:
        if n % 2 == 1:
            r1_1 = line.rstrip()
            filled = True
    if filled:
        filled = False
        count += 1
        # search for motives
        title = r1_0.split(' ')[0]
        seq = r1_1
        pos = search(seq, motive, mismatches)
        if pos > -1:
            fp.write(title + '|POS=' + str(pos) + '\n')
            fp.write(seq + '\n')
            fwd_count += 1
            continue
        else:
            seq = make_revcompl(seq)
            pos = search(seq, motive, mismatches)
            if pos > -1:
                fp.write(title + '|POS=' +str(pos) + '\n')
                fp.write(seq + '\n')
                rev_count += 1
                continue
        no_hit += 1
fp.close()

print("Total seqs: "+str(count)+" Total found: "+str(fwd_count+rev_count)+" (FWD: "+str(fwd_count)+"/REV: "+str(rev_count)+") NO HIT: "+str(no_hit))
