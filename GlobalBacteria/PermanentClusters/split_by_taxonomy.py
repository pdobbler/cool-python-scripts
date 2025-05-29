__author__ = 'avetrot'

import sys
import gzip

fasta_file = sys.argv[1]            # contains md5 in headers
processed_blast = sys.argv[2]
taxonomy_file = sys.argv[3]
category = int(sys.argv[4])
output = sys.argv[5]
treshold = float(sys.argv[6])       # 188.0

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

print("Breakdown to category "+str(category)+" treshol sim+cov "+str(treshold))

# loading of taxonomy...
taxonomy = {}
for line in openfile(taxonomy_file, 'r'):
    parts = line.strip().split('\t')
    name = parts[0]
    categories = parts[1].split(' ')
    if len(categories)>1:
        taxonomy[name] = categories[category].rstrip(';')
print("Taxons loaded - taxonomy: " + str(len(taxonomy)))

# processing blast based on md5
md5_groups = {}
for line in openfile(processed_blast, 'r'):
    parts = line.strip().split('\t')
    if parts[0] == "SEQ_NAME":
        print("header "+line.strip())
    else:
        md5_id = parts[1]
        name = parts[2]
        if name == "NO_HIT":
            md5_groups[md5_id] = "NO_HIT"
        else:
            sim = float(parts[3])
            cov = float(parts[4])
            if taxonomy.has_key(name):
                if sim+cov >= treshold:
                    md5_groups[md5_id] = taxonomy[name]
                else:
                    md5_groups[md5_id] = "NOT_PASS"
            else:
                print("Feature not found: "+name)
print("Blast onfo loaded - groups found: " + str(len(md5_groups)))

filled = False
for n, line in enumerate(openfile(fasta_file)):
    if n % 2 == 0:
        title = line.rstrip()[1:]
    else:
        if n % 2 == 1:
            seq = line.rstrip()
            filled = True
    if filled:
        group_name = md5_groups[title.split('|')[0]]
        with open(output + "_" + group_name + ".fas", 'a') as fp:             # with statement safer than open/close
            fp.write(">" + title + "\n")
            fp.write(seq + "\n")
        filled = False

print("DONE!")
