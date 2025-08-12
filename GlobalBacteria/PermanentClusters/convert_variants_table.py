__author__ = 'vetrot'

import sys
import os
import gzip

variants_file = sys.argv[1] # VARIANTS_TABLE.txt.gz

#################################################

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

##################################################

# 5848bd7129a51af694e2c125db59d9f6        5;8     2;14    -       TACG...
# f65b9f2d4b0d9855c0b01ebe313ef06d  7;19;15 4;5;8 4 TACG...

# CREATE TABLE IF NOT EXISTS `variants` (
#   `id` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL,
#   `hash` varchar(32) NOT NULL,
#   `sequence` TEXT NOT NULL
# );

# CREATE TABLE IF NOT EXISTS `samplevar` (
#   `id` bigint(20) unsigned NOT NULL,
#   `variant` int(10) unsigned NOT NULL,
#   `sample` int(10) unsigned NOT NULL,
#   `abundance` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL
# );

fp_variants = gzip.open("VARIANTS_variants.txt.gz", 'wb')
fp_samplevar = gzip.open("VARIANTS_samplevar.txt.gz", 'wb')
sv_id = 1
var_id = 1
for line in openfile(variants_file):
    vals = line.rstrip().split('\t')
    md5 = vals[0]
    sn = vals[1].split(';') #sample names
    sa = vals[2].split(';') #sample abundances
    cl_id = '0'
    if not vals[3] == '-':
        cl_id = vals[3]
    fp_variants.write(str(var_id) + '\t' + cl_id + '\t' + md5 + '\t' + vals[4] + '\n')
    for x in range(len(sn)):
        fp_samplevar.write(str(sv_id) + '\t' + str(var_id) + '\t' + sn[x] + '\t' + sa[x] + '\t' + cl_id + '\n')
        sv_id += 1
    var_id += 1
fp_samplevar.close()
fp_variants.close()

print("Done :]")
