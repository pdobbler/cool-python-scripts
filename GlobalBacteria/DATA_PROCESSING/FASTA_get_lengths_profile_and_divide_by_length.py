__author__ = 'Wietrack 2025'

import sys
import gzip

fasta_file = sys.argv[1]
cut_off_min = int(sys.argv[2])
cut_off_max = int(sys.argv[3])

################################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
################################################

lenghts_count = {}
loaded = False
titleRead = False

i = 0
n = 0
k = 0

# Open output files in text mode with gzip
first_part = fasta_file[:fasta_file.find('.')]
out_file_under = gzip.open(first_part + "_under_" +str(cut_off_min) + ".fa.gz", 'wt')
out_file_ok = gzip.open(first_part + "_ok_" +str(cut_off_min) + "_" + str(cut_off_max) + ".fa.gz", 'wt')
out_file_above = gzip.open(first_part + "_above_" +str(cut_off_max) + ".fa.gz", 'wt')

#load fasta sequences
for line in openfile(fasta_file):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            l = len(seq)
            # count the lengths
            if lenghts_count.has_key(l):
                lenghts_count[l] = lenghts_count[l] + 1
            else:
                lenghts_count[l] = 1
            # filter
            if l <= cut_off_max and l >= cut_off_min:
                out_file_ok.write(">" + title +'\n')
                out_file_ok.write(seq +'\n')
                k = k + 1
            else:
                if l > cut_off_max:
                    out_file_above.write(">" + title +'\n')
                    out_file_above.write(seq +'\n')
                    i = i + 1
                else:
                    out_file_under.write(">" + title +'\n')
                    out_file_under.write(seq +'\n')
                    n = n + 1

out_file_under.close()
out_file_ok.close()
out_file_above.close()

#save the lengths profile
out_file_profile = open(fasta_file + "_lengths_profile.txt", 'w')
s = set()
for l in lenghts_count:
    s.add(l)

s = sorted(s)
out_file_profile.write("length\tcounts\n")
for l in s:
    out_file_profile.write(str(l) + "\t" + str(lenghts_count[l]) +'\n')
out_file_profile.close()

print(first_part + " " + str(i + n + k) + " sequences processed by length: " + str(k) + " sequence above (>=) " + str(cut_off_max) + " : " + str(i) + "  and uder (>=) " + str(cut_off_min) +  " : " + str(n))

