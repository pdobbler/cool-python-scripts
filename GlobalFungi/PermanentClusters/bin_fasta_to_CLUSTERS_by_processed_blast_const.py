__author__ = 'avetrot'

import sys
import gzip

FASTA = sys.argv[1]
PROCESSED_BLAST = sys.argv[2]
constant = float(sys.argv[3]) # 197.0 or 198.5

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


#############################################

binned={}
for n, line in enumerate(openfile(PROCESSED_BLAST)):
    if n>0:
        vars = line.rstrip().split('\t')
        if vars[3] != '-':
            sum = float(vars[3])+float(vars[4])
            #print("sum "+str(sum))
            if sum>constant:
                binned[vars[0]] = vars[2]
print("Binned sequences were set...")

fp = open(FASTA+".binned", 'w')
fs = open(FASTA+".notbinned", 'w')
filled = False
for n, line in enumerate(openfile(FASTA)):
    if n % 2 == 0:
        r1_0 = line.rstrip()
    else:
        if n % 2 == 1:
            r1_1 = line.rstrip()
            filled = True
    if filled:
        key = r1_0[1:]
        if binned.has_key(key):
            fp.write('>' + binned[key] + '|' + key + '|BINNED_' + str(constant) + '\n')
            fp.write(r1_1 + '\n')
        else:
            fs.write(r1_0 + '\n')
            fs.write(r1_1 + '\n')
        filled = False
fp.close()
fs.close()

print("Done :)")
