__author__ = 'avetrot'

import sys
import gzip

abundance_file = sys.argv[1]
taxonomy_pairs = sys.argv[2]
output = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# loading of taxonomy...
taxonomy = {}
groups = {}
for line in openfile(taxonomy_pairs, 'r'):
    parts = line.strip().split('\t')
    taxonomy[parts[0]] = parts[1]
    groups[parts[1]] = 0
groups["unknown"] = 0
# set list
group_keys = []
for key in groups:
    group_keys.append(key)

print("Taxons loaded - taxonomy: " + str(len(taxonomy)))

# process abundance file
fp = open(output, 'w')
taxonomy_groups = {}
for index, line in enumerate(openfile(abundance_file, 'r')):
    parts = line.strip().split('\t')
    if index == 0:
        # header
        for i, part in enumerate(parts):
            if i > 0:
                if taxonomy.has_key(parts[i]):
                    taxonomy_groups[i] = taxonomy[parts[i]]
                else:
                    taxonomy_groups[i] = "unknown"
    else:
        # reset groups
        for key in group_keys:
            groups[key] = 0
        # get values
        for i, part in enumerate(parts):
            if i > 0:
                groups[taxonomy_groups[i]] += int(parts[i])
        # save
        line = parts[0]
        for key in group_keys:
            line += "\t" + str(groups[key])
        fp.write(line+"\n")
fp.close()
print("DONE!")


