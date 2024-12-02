import sys

table_file = sys.argv[1]
fixed_columns = int(sys.argv[2])    #2    1)ctg_name 2)ctg_length
multi_const = int(sys.argv[3])      #100  for %
out_file = sys.argv[4]


#get col sums....
sums = []
for n, line in enumerate(open(table_file)):
    if n ==0:
        i=0
        vals = line.strip().split("\t")
        for val in vals:
            if i>=fixed_columns:
                #print str(i-fixed_columns)
                sums.append(0)#[i-fixed_columns]=0
            i=i+1
    else:
        i=0
        vals = line.strip().split("\t")
        for val in vals:
            if i>=fixed_columns:
                sums[i-fixed_columns]=sums[i-fixed_columns]+float(val)
            i=i+1

#show sums...
for sum in sums:
    print str(sum)
#normalise table and save...
fp = open(out_file, 'w')
header = ""
for n, line in enumerate(open(table_file)):
    if n ==0:
        fp.write(line.strip()+"\n")
    else:
        vals = line.strip().split("\t")
        new_line = ''
        i=0
        for val in vals:
            if (i>=fixed_columns)and(sums[i-fixed_columns]>0):
                new_line = new_line+str(float(val)/sums[i-fixed_columns]*multi_const)+"\t"
            else:
                new_line = new_line+val+"\t"
            i=i+1
        fp.write(new_line.strip()+"\n")
fp.close()

print "Done :)"