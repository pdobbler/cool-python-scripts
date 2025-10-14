__author__ = 'vetrot'

import sys
import gzip

found_variants = sys.argv[1]  # VARIANTS_variants_GTDB_EXACT.txt
md5_derep_map = sys.argv[2]  # ssu_all_final_for_blast_derep.tab
out_table = sys.argv[3]  # GB1_GTDB_variants.txt

def openfile(filename, mode='r'):
    # For Py3 use 'rt'/'wt' with encoding if needed.
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

####################################################
# id      cl_id   hash                                    seq
###############################################################
# 25937   25335   fe064385e7c0eb5571fe862ba2159af5        TACGTA...
# 50390   4890    e8d2ea834a9fd844d2677f2070a13357        TACGGA...
# 51317   58468   59b845ce8bbc91764f1d7eba585afe0d        TACGTA

md5_vars = {}
for line in openfile(found_variants):
    vals = line.rstrip().split('\t')
    md5_vars[vals[2]] = vals[0]

print("variants loaded")

# 07d81e99fa8897d23f5aa9a86768c692        RS_GCF_000219875.1~NC_017167.1-#6|POS=506|POS=781|POS=506
# 9b0889d6e72e0b4bdc0d8247e54e0876        GB_GCA_007376145.1~VMFP01000020.1|POS=269|POS=543|POS=269
# f7af447c2570b8133a9bf87c403cc53c        GB_GCA_001823385.1~MHUJ01000010.1|POS=483|POS=756|POS=483
# f7af447c2570b8133a9bf87c403cc53c        GB_GCA_001825035.1~MHUO01000007.1|POS=483|POS=756|POS=483
# f7af447c2570b8133a9bf87c403cc53c        GB_GCA_001003905.1~LCQA01000007.1|POS=483|POS=756|POS=483
# f7af447c2570b8133a9bf87c403cc53c        GB_GCA_001003785.1~LCPU01000002.1|POS=483|POS=756|POS=483

vargtdb_id_set = set()
n = 0
fp = open(out_table, 'w')
for line in openfile(md5_derep_map):
    vals = line.rstrip().split('\t')
    if md5_vars.has_key(vals[0]):
        parts = vals[1].split('~')[0].split('_')
        gtdb_id = parts[1] + '_' + parts[2]
        var_id = str(md5_vars[vals[0]])
        vargtdb_id = var_id + '-' +  gtdb_id
        if not vargtdb_id in vargtdb_id_set:
            n += 1
            fp.write(str(n) + '\t' + var_id + '\t' +  gtdb_id + '\n')
            vargtdb_id_set.add(vargtdb_id)
fp.close()

print("Done :)")


