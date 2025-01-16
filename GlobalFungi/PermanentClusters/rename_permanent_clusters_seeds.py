__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip

fasta_file = sys.argv[1]
num_digits = int(sys.argv[2]) #8

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# >S1_Doth_CL164312|9514f5998c5982ee5a0455a0954f55f9|V_2|S_1|P_1|r_2.99116499595e-07|SEED
fp = open(fasta_file+".renamed", 'w')
c = 0
last_rank = float(1000000)
for line in openfile(fasta_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            vals = title.split('|')
            c += 1
            clnamme = "PC"
            for i in range(num_digits-len(str(c))):
                clnamme += '0'
            clnamme += str(c)
            #print(clnamme)
            # save renamed
            fp.write(">" + clnamme + '|' + vals[1] + '|' + vals[2] + '|' + vals[3] + '|' + vals[4] + '|' + vals[5] + "\n")
            fp.write(seq + "\n")
            # count samples sequences
            rank = float(vals[5].split('_')[1])
            if last_rank<rank:
                print("Wrong rank order "+vals[5].split('_')[0]+" last rank "+str(last_rank)+" "+title)
            last_rank = rank
fp.close()


print("Sequences renamed :)")


