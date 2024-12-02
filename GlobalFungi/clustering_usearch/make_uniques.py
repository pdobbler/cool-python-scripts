__author__ = 'avetrot'

import sys
import operator
import gzip

FASTA_in = sys.argv[1]

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################

variants={}
filled = False
for n, line in enumerate(openfile(FASTA_in)):
    if n % 2 == 0:
        r1_0 = line.rstrip()
    else:
        r1_1 = line.rstrip()
        filled = True
    if filled:
        if variants.has_key(r1_1):
            variants[r1_1] += 1
        else:
            variants[r1_1] = 1
        filled = False
print("Variants counts were set...")

# >Uniq1;size=320;
#x = {1: 2, 3: 4, 4: 3, 2: 1, 0: 0}
sorted_variants = sorted(variants.items(), key=operator.itemgetter(1), reverse=True)
print("Variants sorted...")

fp = open(FASTA_in+".uniq", 'w')
index = 1
for v in sorted_variants:
    #print('size '+str(variants[v]))
    #print('seq  ' + v[0])
    #print('size ' + str(v[1]))
    fp.write('>Uniq' + str(index) + ';size=' + str(v[1]) + '\n')
    fp.write(v[0] + '\n')
    index += 1
fp.close()

print("Done :)")
