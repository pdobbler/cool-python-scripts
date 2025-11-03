__author__ = 'avetrot'

import sys
import gzip

clusters_tax = sys.argv[1] # TAXONOMY_CLUSTERS.sorted.txt
variants_tab = sys.argv[2] # VARIANTS_variants.txt.gz
samplevars_tab = sys.argv[3] # VARIANTS_samplevar.txt.gz
samples_tab = sys.argv[4] # SAMPLES_BASIC3.txt
out_tab = sys.argv[5]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# 1       GB00000001.1    Bradyrhizobium sp000244915      Bradyrhizobium_503372   100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rhizobiales_505101; f__Xanthobacteraceae; g__Bradyrhizobium_503372; s__Bradyrhizobium sp000244915    6472eb8b1e09f892aca2f23182962903
cl_seeds = {}
for line in openfile(clusters_tax):
    parts = line.strip().split('\t')
    cl_seeds[parts[7]] = -1 # 6472eb8b1e09f892aca2f23182962903 vs 1

print("Clusters taxonomy loaded: " + str(len(cl_seeds)))

# zcat VARIANTS_variants.txt.gz | head
# CREATE TABLE IF NOT EXISTS `variants` (
#   `id` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL,
#   `hash` varchar(32) NOT NULL,
#   `sequence` TEXT NOT NULL
# );
# 1       7933    7b26301f5893e689b354a0a981910c02        TACGGAGGGTGCAATCGTTATTCGGATTTACTGGGCTTAAAGCGCGCGTAGGCGGCCGTTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCCTGGAAGTGCATTTGATACTATTCGGCTTGAGTATGGGAGAGGGAAGTGGAATTCCTGGTGTAGAGGTGAAATTCGTAGATATCAGGAGGAACACCGGTGGCGAAGGCGACTTCCTGGACCAATACTGACGCTGAGGCGCGAAGGCGTGGGGAGCAAACAGG
# 2       375     99a556f970834de4c7d9eda6d9a4b423        GACGAACCGTGCGAACGTTGTTCGGAATCACTGGGCTTAAAGGGCGCGTAAGCGGCTTGCCAAGTCAGTGGTGAAATCCCGCAGCTTAACTGCGGAAGTGCCTTTGATACTGGCGAGCTCGAGGGAGGTAGGGGTATGTGGAACTTCTGGTGGAGCGGTGAAATGCGTTGATATCAGAAGGAACGCCGGTGGCGAAAGCGACGTACTGGACCTCTTCTGGCGCTGAGGCGCGAAAGCTAGGGGAGCAAACGGG
vars_seq_count = {}
vars_samples = {}
vars_papers = {}
cl_vars = {}
for line in openfile(variants_tab):
    parts = line.strip().split('\t')
    md5 = parts[2]
    cl_id = parts[1]
    var_id = parts[0]
    # variants in each cluster
    if cl_vars.has_key(cl_id):
        cl_vars[cl_id] += 1
    else:
        cl_vars[cl_id] = 1
    # seed var 
    if cl_seeds.has_key(md5):
        vars_seq_count[var_id] = 0
        vars_samples[var_id] = set()
        vars_papers[var_id] = set()
        cl_seeds[md5] = var_id

# head SAMPLES_BASIC3.txt
# CREATE TABLE IF NOT EXISTS `samples_basic` (
#   `id` int NOT NULL PRIMARY KEY,
#   `paper` int NOT NULL,
# 3253    98      GB01023803S     soil    -40.33  175.24  Australia       2023    2023    forest  515F/806R       13.9    990.0   5.3     11.78   70620   0
# 3382    88      GB01020485S     soil    52.377418       14.290489       Europe  2023    2023    cropland        341F/785R       9.3     546.0   6.0     5.12    88204   0
# 10120   72      GB01016029S     soil    -54.902614      -68.011204      South America   2016    2016    wetland 515F-Y/928R     5.6     589.0   5.6     16.53   36897   0
sample_papers = {}
for line in openfile(samples_tab):
    parts = line.strip().split('\t')
    sample_papers[parts[0]] = parts[1]

# zcat VARIANTS_samplevar.txt.gz | head
# CREATE TABLE IF NOT EXISTS `samplevar` (
#   `id` bigint(20) unsigned NOT NULL,
#   `variant` int(10) unsigned NOT NULL,
#   `sample` int(10) unsigned NOT NULL,
#   `abundance` int(10) unsigned NOT NULL,
#   `cl_id` int(10) unsigned NOT NULL
# );
# 1       1       6854    1       7933
# 2       2       1073    1       375
cl_seq_count = {}
cl_samples = {}
cl_papers = {}
for line in openfile(samplevars_tab):
    parts = line.strip().split('\t')
    cl_id = parts[4]
    sample = parts[2]
    var = parts[1]
    seqs = int(parts[3])
    # samples counts
    if cl_samples.has_key(cl_id):
        samples = cl_samples[cl_id]
        samples.add(sample)
        cl_samples[cl_id] = samples
        # papers
        papers = cl_papers[cl_id]
        papers.add(sample_papers[sample])
        cl_papers[cl_id] = papers
    else:
        samples = set()
        samples.add(sample)
        cl_samples[cl_id] = samples
        # papers
        papers = set()
        papers.add(sample_papers[sample])
        cl_papers[cl_id] = papers
    # cluster seq counts
    if cl_seq_count.has_key(cl_id):
        cl_seq_count[cl_id] += seqs        
    else:
        cl_seq_count[cl_id] = seqs
    # seed seq counts
    if vars_seq_count.has_key(var):
        vars_seq_count[var] += seqs
        # seed seq samples
        vars_samples[var].add(sample)
        vars_papers[var].add(sample_papers[sample])


# write output
fp = open(out_tab, 'w')
line = "cl_id\tcl_name\tps_name\tgen_name\tsim\tcov\ttaxonomy\tseed_md5"
line += "\tcl_samples\tcl_papers\tcl_seqs"
line += "\tseed_samples\tseed_papers\tseed_seqs"
line += "\tcl_variants"
fp.write(line + '\n')
for line in openfile(clusters_tax):
    parts = line.strip().split('\t')
    line = line.strip()
    cl_id = parts[0]
    md5 = parts[7]
    var = cl_seeds[md5]
    # cl sample counts and sequences
    line += "\t" + str(len(cl_samples[cl_id])) + "\t" + str(len(cl_papers[cl_id])) + "\t" + str(cl_seq_count[cl_id])
    # var sample counts and sequences
    line += "\t" + str(len(vars_samples[var])) + "\t" + str(len(vars_papers[var])) + "\t" + str(vars_seq_count[var])
    # cluster var counts
    line += "\t" + str(cl_vars[cl_id])
    fp.write(line + '\n')
fp.close()

print("DONE :)")
    




