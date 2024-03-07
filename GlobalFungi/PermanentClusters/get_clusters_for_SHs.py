__author__ = 'avetrot'

import sys

SH_in_CLUSTERS = sys.argv[1]

shs_clusters = {}
for line in open(SH_in_CLUSTERS):
    values = line.rstrip().split('\t')
    cl_name = values[0]
    shs = values[1].split(';')
    for sh in shs:
        if shs_clusters.has_key(sh):
            clusters = shs_clusters[sh]
            clusters.add(cl_name)
            shs_clusters[sh] = clusters
        else:
            clusters = set()
            clusters.add(cl_name)
            shs_clusters[sh] = clusters

print("file loaded...")

fp = open("clusters_of_shs.txt", 'w')
for sh in shs_clusters:
    clusters = shs_clusters[sh]
    fp.write(sh + '\t' + ';'.join(clusters) + '\t' + str(len(clusters)) + '\n')
fp.close()

print('DONE :]')

