__author__ = 'vetrot'

import sys
from random import randint
import array

size = int(sys.argv[1])
fasta_file = sys.argv[2]
subname = sys.argv[3]

open(subname, "w")

filled = False

num_seqs = sum(1 for line in open(fasta_file))/2
print "number of fasta sequences "+str(num_seqs)
a = array.array('i',(i for i in range(0,num_seqs)))
print "processing random choose from array of lenght "+str(len(a))
random_list = {}
for i in range(0,size):
    if len(a)>0:
        index = randint(0, len(a)-1)
        if i % 10000 == 0:
            print "random progress is "+str(i)
        # print "index to remove "+str(index)
        random_list[a[index]] = 0
        del a[index]
print "random chooose done..."+str(len(random_list))

# check tag and save to new R1 and R2 file
def save_fasta(r1_0,r1_1):
    with open(subname, "a") as fileOUT:
        # print r1_0 + " intact", best_pos
        fileOUT.write('%s\n' % r1_0)
        fileOUT.write('%s\n' % r1_1)
    return;

index = 0
for n, line1 in enumerate(open(fasta_file)):
    if n % 200000 == 0:
        print n / 2
    if n % 2 == 0:
        r1_0 = line1.rstrip()
    else:
        if n % 2 == 1:
            r1_1 = line1.rstrip()
            filled = True
    if filled:
        if random_list.has_key(index):
            save_fasta(r1_0,r1_1)
        index = index + 1
        filled = False

print "Done."