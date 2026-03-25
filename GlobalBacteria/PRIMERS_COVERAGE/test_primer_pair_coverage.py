__author__ = 'avetrot'

import sys
import gzip

FASTA = sys.argv[1]
pair_name = sys.argv[2]
fwd_primer_seq = sys.argv[3]
rev_primer_seq = sys.argv[4]

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
    }.get(x, 'N')    # 'N' is default if x not found

def make_revcompl(seq):
    rev_compl_seq = ''
    for i in range(len(seq)):
        rev_compl_seq = rev_compl_seq + compl_nucl(seq[len(seq) - (i + 1)].upper())
    return rev_compl_seq

def compare_strings_limit(s1, s2, max_mismatch):
    mm = 0
    for a, b in zip(s1, s2):
        if a != b:
            mm += 1
            if mm > max_mismatch:
                return mm
    return mm


def search(sequence, adapt, mismatch):
    n = len(sequence)
    m = len(adapt)

    if m == 0:
        return 0
    if m > n:
        return -1

    # Fast path for exact match
    if mismatch == 0:
        return sequence.find(adapt)

    limit = n - m + 1
    for i in range(limit):
        mm = 0
        for j in range(m):
            if sequence[i + j] != adapt[j]:
                mm += 1
                if mm > mismatch:
                    break
        if mm <= mismatch and m > mismatch + 1:
            return i

    return -1

DEGENERATE_MAP = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

########################
# OBJECT PRIMER
########################

class PrimerVariants(object):

    def __init__(self, fwd_primer, rev_primer, name):
        self.fwd_primer = fwd_primer.upper()
        self.rev_primer = rev_primer.upper()
        self.fwd_variants = self._generate_variants(self.fwd_primer)
        self.rev_variants = self._generate_variants(self.rev_primer)
        self.fwd_hits = dict((fwd_variant, []) for fwd_variant in self.fwd_variants)
        self.rev_hits = dict((rev_variant, []) for rev_variant in self.rev_variants)
        self.taxa = {}
        self.name = name

    def _generate_variants(self, primer):
        variants = ['']
        for base in primer:
            if base not in DEGENERATE_MAP:
                raise ValueError("Unknown base in primer: %s" % base)

            new_variants = []
            for prefix in variants:
                for nucleotide in DEGENERATE_MAP[base]:
                    new_variants.append(prefix + nucleotide)
            variants = new_variants
        return variants

    def reset(self):
        self.fwd_hits = dict((fwd_variant, 0) for fwd_variant in self.fwd_variants)
        self.rev_hits = dict((rev_variant, 0) for rev_variant in self.rev_variants)
        self.taxa = {}

    def count_hits_in_fasta(self, fasta_file):
        self.reset()
        filled = False
        for n, line in enumerate(openfile(fasta_file)):
            if n % 2 == 0:
                r1_0 = line[1:].rstrip()
            else:
                if n % 2 == 1:
                    r1_1 = line.rstrip()
                    filled = True
            if filled:
                filled = False
                taxa = r1_0.split(';s__')[0]
                seq = r1_1
                found = -1

                for fwd_variant in self.fwd_variants:
                    pos = seq.find(fwd_variant)
                    if pos > -1:
                        self.fwd_hits[fwd_variant] += 1
                        found = 1
                        break

                if found == -1:
                    for fwd_variant in self.fwd_variants:
                        pos = seq.find(make_revcompl(fwd_variant))
                        if pos > -1:
                            self.fwd_hits[fwd_variant] += 1
                            found = 2
                            break

                # check reverse
                if found == 2:
                    for rev_variant in self.rev_variants:
                        pos = seq.find(rev_variant)
                        if pos > -1:
                            self.rev_hits[rev_variant] += 1
                            found = 3
                            break

                if found == 1:
                    for rev_variant in self.rev_variants:
                        pos = seq.find(make_revcompl(rev_variant))
                        if pos > -1:
                            self.rev_hits[rev_variant] += 1
                            found = 4
                            break

                if found > 2:
                    if self.taxa.has_key(taxa):
                        self.taxa[taxa] += 1
                    else:
                        self.taxa[taxa] = 1

    def get_fwd_hits(self):
        return self.fwd_hits

    def get_rev_hits(self):
        return self.rev_hits

    def get_fwd_variants(self):
        return self.fwd_variants

    def get_rev_variants(self):
        return self.rev_variants

    def save_taxa_coverage(self):
        safe_name = self.name.replace('/', '_')
        fp = open(safe_name + "_taxonomy_coverage.txt", 'w')
        for t in self.taxa:
            fp.write(t + "\t" + str(self.taxa[t]) + "\n")
        fp.close()

    def __repr__(self):
        return "PrimerVariants(primer='%s', variants=%d)" % (
            self.primer, len(self.variants)
        )

########################
# EXECUTION
########################

pv = PrimerVariants(fwd_primer_seq, rev_primer_seq, pair_name)
print("FWD")
print pv.get_fwd_variants()
print("REV")
print pv.get_rev_variants()
pv.count_hits_in_fasta(FASTA)
pv.save_taxa_coverage()
print("FWD HITS")
print pv.get_fwd_hits()
print("REV HITS")
print pv.get_rev_hits()

print("Done")


