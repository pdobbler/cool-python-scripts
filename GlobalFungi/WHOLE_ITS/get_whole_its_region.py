__author__ = 'vetrot'

import sys

derep_fasta = sys.argv[1]        # test_dereplicated.fa
itsx_positions = sys.argv[2]     # itsx_test_dereplicated.positions.txt
complete_ITS_fasta = sys.argv[3]

# seq_name  630 bp.   SSU: 1-80   ITS1: 81-630   5.8S: Not found   ITS2: Not found   LSU: Not found   Broken or partial sequence, no 5.8S!
# g0010484|size=1	577 bp.	SSU: 1-80	ITS1: 81-242	5.8S: 243-399	ITS2: 400-545	LSU: 546-577
ranges = {}
done = {}
i = 0
fp = open("suspicious_positions.txt", 'w')
for line in open(itsx_positions):
    parts = line.rstrip().split('\t')
    seq_name = parts[0]
    #
    SSU = parts[2].split(': ')[1]
    ITS1 = parts[3].split(': ')[1]
    S58 = parts[4].split(': ')[1]
    ITS2 = parts[5].split(': ')[1]
    LSU = parts[6].split(': ')[1]
    #
    if SSU != "Not found" and ITS1 != "Not found" and S58 != "Not found" and ITS2 != "Not found" and LSU != "Not found":
        ranges[seq_name] = ITS1.split('-')[0] + "-" + ITS2.split('-')[1]
        done[seq_name] = line.rstrip()
        if ITS1 == "No end" or S58 == "No end" or ITS2 == "No end":
            i += 1
            fp.write(line.rstrip() + '\n')
fp.close()
print("Position for complete region were loaded... ("+str(len(ranges))+") suspicious positions (" + str(i) + ")")
if i > 0:
    print("For suspicious positions - check suspicious_positions.txt")

i = 0
filled = False
fp = open(complete_ITS_fasta, 'w')
for n, line in enumerate(open(derep_fasta)):
    if n % 2 == 0:
        title = line[1:].rstrip()
    else:
        sequence = line.rstrip()
        filled = True

    if filled:
        filled = False
        if ranges.has_key(title):
            done[title] = ""
            i += 1
            pos = ranges[title].split('-')
            fp.write('>' + title + '\n')
            fp.write(sequence[int(pos[0])-1:int(pos[1])] + '\n')
fp.close()

print("Complete sequences were extracted " + str(i) + " out of "+ str(len(ranges)))
if i != len(ranges):
    print("Error - some sequences are missing!!! (" + str(len(ranges)-i)+") - check not_found_positions.txt")
else:
    print("Extraction was successful!")

fp = open("not_found_positions.txt", 'w')
for title in done:
    if done[title] != "":
        fp.write(done[title] + '\n')
fp.close()

print("DONE:]")