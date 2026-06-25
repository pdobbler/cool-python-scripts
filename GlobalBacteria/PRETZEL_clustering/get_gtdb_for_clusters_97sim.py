__author__ = 'vetrot'

import sys
import gzip

cluster_vs_gtdb16S = sys.argv[1]  # RAW_GTDB_CLUSTERS97sim.txt
out_table = sys.argv[2]  # GB2_GTDB_CLUSTERS97sim.txt

def openfile(filename, mode='r'):
    # For Py3 use 'rt'/'wt' with encoding if needed.
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def hit_numeric_key(hit):
    """
    Extract the main numeric ID from a hit name like 'GB00002469.1' -> 2469.
    If multiple numeric groups are present, ignore the final version suffix.
    Return '-' if hit == 'NO_HIT'.
    """
    import re
    if hit == "NO_HIT":
        return '-'

    # Match GB + digits + .version
    m = re.match(r'^[A-Z]+0*(\d+)\.\d+$', hit)
    if m:
        return int(m.group(1))

    # Fallback: take first integer group
    m = re.search(r'(\d+)', hit)
    return int(m.group(1)) if m else 10**12

# CL00000001      GB_GCA_000013685.1
# CL00000001      GB_GCA_000219645.1
# CL00000001      GB_GCA_000374205.1
# CL00000001      GB_GCA_000465325.1
# CL00000001      GB_GCA_000472765.1
# CL00000001      GB_GCA_000473005.1
# CL00000001      GB_GCA_001549695.1

fp = open(out_table, 'w')
n = 0
for line in openfile(cluster_vs_gtdb16S):
        vals = line.rstrip().split('\t')
        cluster_id = hit_numeric_key(vals[0])
        gtdb_id = vals[1]
        n += 1
        fp.write(str(n) + '\t' + str(cluster_id) + '\t' +  gtdb_id + '\n')
fp.close()

print("Done :)")

