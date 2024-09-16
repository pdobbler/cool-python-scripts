__author__ = 'vetrot'

import sys
import ntpath
import gzip

in_file = sys.argv[1]
group_size = int(sys.argv[2])

#############################################
# GZIP OPENING
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#############################################

part = 1
k = 1

file_name = ntpath.basename(in_file).rstrip().split('.')[0]
fp = open(file_name+str(part)+'.fas', 'w')

for n, line in enumerate(openfile(in_file)):
    ch = line[0]
    if ch == '>':
        if ((k % group_size)==0):
            part = part+1
            fp.close()
            print "Part: "+str(part)
            fp = open(file_name+str(part)+'.fas', 'w')
        k=k+1
    fp.write(line)
fp.close()
print "Done."

