__author__ = 'avetrot'

import sys
import gzip

processed_blast = sys.argv[1]
taxonomy_file = sys.argv[2]
category = int(sys.argv[3])
output = sys.argv[4]
treshold = float(sys.argv[5])

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
        taxonomy[name] = categories[category]
print("Taxons loaded - taxonomy: " + str(len(taxonomy)))

# processing blast
tax_counts = {}
for line in openfile(processed_blast, 'r'):
    parts = line.strip().split('\t')
    if parts[0] == "SEQ_NAME":
        print("header "+line.strip())
    else:
        name = parts[2]
        if name == "NO_HIT":
            if tax_counts.has_key("NO_HIT"):
                tax_counts["NO_HIT"] += 1
            else:
                tax_counts["NO_HIT"] = 1
        else:
            sim = float(parts[3])
            cov = float(parts[4])
            if taxonomy.has_key(name):
                if sim+cov >= treshold:
                    taxon = taxonomy[name]
                    if tax_counts.has_key(taxon):
                        tax_counts[taxon] += 1
                    else:
                        tax_counts[taxon] = 1
                else:
                    if tax_counts.has_key("NOT_PASS"):
                        tax_counts["NOT_PASS"] += 1
                    else:
                        tax_counts["NOT_PASS"] = 1
            else:
                print("Feature not found: "+name)
print("Blast onfo loaded - tax_counts: " + str(len(tax_counts)))

fp = open(output, 'w')
fp.write("Taxon\tcounts\n")
for taxon in tax_counts:
    fp.write(taxon+"\t"+str(tax_counts[taxon])+"\n")
fp.close()
print("DONE!")
