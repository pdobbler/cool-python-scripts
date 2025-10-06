__author__ = 'vetrot'

import sys
import operator
import collections
import gzip
import hashlib

samplevar_file = sys.argv[1] # VARIANTS_samplevar.txt.gz
sample_pairs = sys.argv[2]   # VARIANTS_TABLE_SAMPLE_PAIRS.txt
cl_max = int(sys.argv[3])    # 100000

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def make_cluster_name(num):
    """
    Generate cluster name in the format 'GB00000001.1'
    from a given integer number.
    Works in Python 2.7.
    """
    return "GB{:08d}.1".format(num)

# SampleName      SampleID        occurence
# 16349   GB01007340S     5521
# 9601    GB01012613S     19312

samples = {}
samples_occ = {}
for n, line in enumerate(openfile(sample_pairs)):
    if n > 0:
        vals = line.rstrip().split('\t')
        samples[vals[0]] = vals[1]
        samples_occ[vals[1]] = float(vals[2])


# id      var_id  sam_id  sam_ab  cl_id
# 18      13      96      1       195
# 19      13      9727    1       195
# 20      13      9214    1       195
# 21      14      2076    1       1000000000000
# 22      15      454     1       281
# 23      16      1167    1       6346
# fp_samplevar.write(str(sv_id) + '\t' + str(var_id) + '\t' + sn[x] + '\t' + sa[x] + '\t' + cl_id + '\n')

cl_samples = {}
sample_names = set()
for line in openfile(samplevar_file):
    vals = line.rstrip().split('\t')
    cl_id = int(vals[4])
    if cl_id <= cl_max:
        sam_name = samples[vals[2]]
        sam_ab = int(vals[3])
        sample_names.add(sam_name)
        if cl_samples.has_key(cl_id):
            sample_ab = cl_samples[cl_id]
            if sample_ab.has_key(sam_name):
                sample_ab[sam_name] += sam_ab
            else:
                sample_ab[sam_name] = sam_ab
            cl_samples[cl_id] = sample_ab
        else:
            sample_ab = {}
            sample_ab[sam_name] = sam_ab
            cl_samples[cl_id] = sample_ab

print("clusters: "+str(len(cl_samples)))

sample_list = []
for sam_name in sample_names:
    sample_list.append(sam_name)

print("sample names: "+str(len(sample_list)))

fp = open("OTUTABLE_TOP_CLUSTERS_"+str(cl_max)+".txt", 'w')
line = "GB_CLUSTER"
for sam_name in sample_list:
    line += sam_name
fp.write(line +'\n')
for i in range(1, cl_max + 1):
    sample_ab = cl_samples[i]
    line = make_cluster_name(i)
    for sam_name in sample_list:
        if sample_ab.has_key(sam_name):
            line += "\t" + str(float(sample_ab[sam_name])/samples_occ[sam_name]*100)
        else:
            line += "\t0"
    fp.write(line +'\n')
fp.close()

print("Done :)")



