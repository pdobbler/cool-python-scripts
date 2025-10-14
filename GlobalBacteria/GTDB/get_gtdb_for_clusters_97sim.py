__author__ = 'vetrot'

import sys
import gzip

cluster_vs_gtdb16S = sys.argv[1]  # TABLE_3_97sim.txt
out_table = sys.argv[2]  # GB1_GTDB_CLUSTERS97sim.txt

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

################################################################
#GB_md5  count_of_same_MAG_16S_contig_marker     list_of_contigs
#GB00000084.1    6       GB_GCA_035307685.1~DATGHV010000039.1;GB_GCA_035713825.1~DASTMX010000037.1;GB_GCA_035538875.1~DATLJP010000015.1;GB_GCA_035568175.1~DATMON010000009.1;RS_GCF_009394175.1~NZ_RYCI01000193.1;GB_GCA_036825615.1~DAISXY010000006.1
#GB00022730.1    2       GB_GCA_019637225.1~JAHBWC010000011.1;GB_GCA_021462865.1~JAJTYX010000092.1
#GB00458257.1    1       GB_GCA_039795895.1~JBDKKW010000026.1
#GB00132115.1    1       GB_GCA_016713305.1~JADJPK010000007.1

fp = open(out_table, 'w')
i = 0
for line in openfile(cluster_vs_gtdb16S):
    if i>0:
        vals = line.rstrip().split('\t')
        cluster_id = hit_numeric_key(vals[0])
        for gtdb_acc in vals[2].split(';'):
            parts = gtdb_acc.split('~')[0].split('_')
            gtdb_id = parts[1] + '_' + parts[2]
            fp.write(str(i) + '\t' + str(cluster_id) + '\t' +  gtdb_id + '\n')
    i += 1
fp.close()

print("Done :)")