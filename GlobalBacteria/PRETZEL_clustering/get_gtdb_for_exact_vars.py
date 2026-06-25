__author__ = 'vetrot'

import sys
import gzip

found_variants = sys.argv[1]  # VARIANTS_variants_GTDBonly.txt
md5_derep_map = sys.argv[2]  # ssu_all_r232_GB_EXTRACTED_md5_uniq_acc.txt
out_table = sys.argv[3]  # GB2_GTDB_variants.txt

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

# head VARIANTS_variants_GTDBonly.txt
# 12956   5086    c9a2743971d12764732ddc336e25cac8        TAC...
# 26596   38123   a86714dc74c2023ee01017a49444a659        TACG...
# 36056   42747   fe064385e7c0eb5571fe862ba2159af5        TACGTAGG....
# 69885   1612    e8d2ea834a9fd844d2677f2070a13357        TACG...


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

# head ssu_all_r232_GB_EXTRACTED_md5_uniq_acc.txt
# 000033a7eb4277829abd276b0e624f45        GB_GCA_902805585.1
# 0000672cf37091ed24be6595ac42659c        GB_GCA_016212795.1
# 0000d1a50f0c9ae122d1b45df38c4638        GB_GCA_002420765.1
# 0000f9a9a96f81bcfc9ab4393bf33e0e        GB_GCA_016778845.1
# 0000f9a9a96f81bcfc9ab4393bf33e0e        GB_GCA_038141425.1
# 0001d123420b59585627edf5a1292ae8        RS_GCF_001017175.1


vargtdb_id_set = set()
n = 0
fp = open(out_table, 'w')
for line in openfile(md5_derep_map):
    vals = line.rstrip().split('\t')
    if md5_vars.has_key(vals[0]):
        parts = vals[1].split('_')
        gtdb_id = parts[1] + '_' + parts[2]
        var_id = str(md5_vars[vals[0]])
        vargtdb_id = var_id + '-' +  gtdb_id
        if not vargtdb_id in vargtdb_id_set:
            n += 1
            fp.write(str(n) + '\t' + var_id + '\t' +  gtdb_id + '\n')
            vargtdb_id_set.add(vargtdb_id)
fp.close()

print("Done :)")



