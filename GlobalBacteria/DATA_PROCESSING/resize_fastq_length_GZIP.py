__author__ = 'vetrot'

import sys
import gzip

FASTQ_file = sys.argv[1]
OUT_file = sys.argv[2]
length = int(sys.argv[3])


def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

r1_0 = ''
r1_1 = ''
r1_2 = ''
r1_3 = ''
filled = False
fp = gzip.open(OUT_file+'.gz', 'wt')

def save_by_tag(r1_0,r1_1,r1_2,r1_3,length):
    max_len = len(r1_1)
    if length > max_len:
        length = max_len
    fp.write('%s\n' % r1_0)
    fp.write('%s\n' % r1_1[:length])
    fp.write('%s\n' % r1_2)
    fp.write('%s\n' % r1_3[:length])
    return;

for n, line in enumerate(openfile(FASTQ_file)):
    if n % 40000 == 0:
        print n / 4
    if n % 4 == 0:
        r1_0 = line.rstrip()
        #print "line1 %s" % line1
        #print "line2 %s" % line2
    else:
        if n % 4 == 1:
            r1_1 = line.rstrip()
        if n % 4 == 2:
            r1_2 = line.rstrip()
        if n % 4 == 3:
            r1_3 = line.rstrip()
            filled = True

    if filled:
        save_by_tag(r1_0,r1_1,r1_2,r1_3,length)
        filled = False

fp.close()

print "Done."
