__author__ = 'avetrot'

import sys
import gzip
import hashlib

CLUSTERED_VARS = sys.argv[1] #REL4_ITS2_FINAL_qualified_clustered_and_binned.fa
RAW_GF_FASTA = sys.argv[2]   #REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz
sample_list = sys.argv[3]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# GF01000754S

sample_set = set()
for line in openfile(sample_list, 'r'):
    sample_set.add(line.strip())

print("Sample list looaded... "+str(len(sample_set)))


# >PCL000001|ad765e6d39aa1c7a3811388dd6c850c6|V_33910054|S_12224|P_146|r_0.728005994798|100.0 

vars = {}
for line in openfile(CLUSTERED_VARS, 'r'):
    ch = line[0]
    if ch == '>':
        parts = line[1:].strip().split('|')
        vars[parts[1]] = parts[0]

print("Sequence vars loaded... "+str(len(vars)))

# >GF4S07047b|Sun_2021_PK|ERR4887914.1283241

samples_map = {}
papers_map = {}
cl_names = set()
cl_size = 0
total_size = 0
for line in openfile(RAW_GF_FASTA, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            total_size += 1
            hash = hashlib.md5(line.strip().encode()).hexdigest()
            if vars.has_key(hash):
                cl_size += 1
                sample_name = title.split('|')[0]
                paper_name = title.split('|')[1]
                if sample_name in sample_set:
                    cl_name = vars[hash]
                    cl_names.add(cl_name)
                    if cl_name in samples_map:
                        # samples
                        samples = samples_map[cl_name]
                        samples.add(sample_name)
                        samples_map[cl_name] = samples
                        # papers
                        papers = papers_map[cl_name]
                        papers.add(paper_name)
                        papers_map[cl_name] = papers
                    else:
                        # samples
                        samples = set()
                        samples.add(sample_name)
                        samples_map[cl_name] = samples
                        # papers
                        papers = set()
                        papers.add(paper_name)
                        papers_map[cl_name] = papers

print("Sequences processed... total: "+str(total_size)+" clustered: "+str(cl_size)+" non-clustered: "+str(total_size-cl_size))

fp = open(CLUSTERED_VARS+".SAMPLES_and_PAPERS.txt", 'w')
for cl_name in cl_names:
    fp.write(cl_name + '\t' + str(len(samples_map[cl_name])) + '\t' + str(len(papers_map[cl_name])) + '\n')
fp.close()

print("DONE. "+str(len(vars)))

