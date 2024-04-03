__author__ = 'daniel_pusina_morais'

import sys

barcodes = sys.argv[1]
R1 = sys.argv[2]
R2 = sys.argv[3]

# loading of barcodes...
samples = {}
barcode = {}
min = 1000
max = 0
for line in open(barcodes):
    dat = line.rstrip().split('\t')
    samples[dat[2]+'|'+dat[4]] = dat[0]
    barcode[dat[1]+'|'+dat[3]] = dat[2]+'|'+dat[4]
    # set min and max length...
    if len(dat[1])< min:
        min = len(dat[1])
    if len(dat[1])> max:
        max = len(dat[1])
    if len(dat[3])< min:
        min = len(dat[3])
    if len(dat[3])> max:
        max = len(dat[3])

print("min "+ str(min))
print("max "+ str(max))

# reading of fastq R1...

#with open(R1) as f1, open(R2) as f2:
#    while True:
#        line1 = f1.readline()
#        line2 = f2.readline()
#        if line1 == '':
#            break
#        print(line1)
#        print(line2)



filled = False
n = 0
fp1 = open(R1 + ".cut", 'w')
fp2 = open(R2 + ".cut", 'w')
with open(R1) as f1, open(R2) as f2:
    while True:
        line1 = f1.readline().rstrip()
        line2 = f2.readline().rstrip()
        if line1 == '':
            break
        if n % 400000 == 0:
            print n / 4
        if n % 4 == 0:
            r1_0 = line1
            r2_0 = line2
        else:
            if n % 4 == 1:
                r1_1 = line1 # nucleotides 1
                r2_1 = line2  # nucleotides 2
            if n % 4 == 2:
                r1_2 = line1
                r2_2 = line2
            if n % 4 == 3:
                r1_3 = line1 # quality 1
                r2_3 = line2  # quality 2
                filled = True
        if filled:
            #i = test_fasta(r1_0,r1_1,r1_2,r1_3,i)
            sample = ""
            for i in range(min, max+1):
                bar1 = r1_1[0:i]
                for j in range(min, max+1):
                    bar2 = r2_1[0:j]
                    found = False
                    #normal = False
                    # testing in normal order...
                    if barcode.has_key(bar1+'|'+bar2):
                        sample = barcode[bar1+'|'+bar2]
                        #print(sample + " normal -> " + samples[sample])
                        found = True
                        normal = True
                    # reverse order...
                    if barcode.has_key(bar2+'|'+bar1):
                        sample = barcode[bar2+'|'+bar1]
                        #print(sample + " reverse -> " + samples[sample])
                        found = True
                        normal = False
                    if found:
                        r1_0 = r1_0.split(' ')[0]
                        r2_0 = r2_0.split(' ')[0]
                        if normal:
                            # cut the barcode and save fastq...
                            fp1.write(r1_0 + '|' + samples[sample] + ' 1:MOD' + '\n')
                            fp1.write(r1_1[i:] + '\n')
                            fp1.write(r1_2 + '\n')
                            fp1.write(r1_3[i:] + '\n')
                            fp2.write(r2_0 + '|' + samples[sample] + ' 2:MOD' + '\n')
                            fp2.write(r2_1[j:] + '\n')
                            fp2.write(r2_2 + '\n')
                            fp2.write(r2_3[j:] + '\n')
                        else:
                            # cut the barcode and save fastq in correct order...
                            fp2.write(r1_0 + '|' + samples[sample] + ' 2:MOD' + '\n')
                            fp2.write(r1_1[i:] + '\n')
                            fp2.write(r1_2 + '\n')
                            fp2.write(r1_3[i:] + '\n')
                            fp1.write(r2_0 + '|' + samples[sample] + ' 1:MOD' + '\n')
                            fp1.write(r2_1[j:] + '\n')
                            fp1.write(r2_2 + '\n')
                            fp1.write(r2_3[j:] + '\n')
                        # Break the inner loop...
                        break
                else:
                    # Continue if the inner loop wasn't broken.
                    continue
                # Inner loop was broken, break the outer.
                break
            # ready for another sequence...
            filled = False
        n += 1

fp1.close()
fp2.close()

#cmd = 'usearch -cluster_otus '+fasta_file+'.uniques.fa -minsize 1 -otus '+fasta_file+'.otus.fa -relabel Otu -uparseout '+fasta_file+'.uparse.txt'
#print(cmd)
#stream = os.popen(cmd)
#output = stream.read()
#print(output)


print "done :)"
