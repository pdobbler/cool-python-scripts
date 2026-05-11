__author__ = 'avetrot'

import sys
import gzip

cl_vars_file = sys.argv[1]  # GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz
processed_blast = sys.argv[2] # CLUSTERS_IDENT_greenegenes2_2024_09_PROCESSED.txt
all_seqs_fasta = sys.argv[3]  # GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz
samplesAndStud = sys.argv[4] # samples_and_studies.txt
output = sys.argv[5]

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

#######################
# 1. CLUSTERS SAMPLES #
#######################

# >CL00000001|d829bee4984f82ffc2453212157caf96;samples=25791;relabund_sum=77.2086052524;size=5579362|100.0
# >CL00000001|6472eb8b1e09f892aca2f23182962903;samples=23821;relabund_sum=111.2113833946;size=6712798|99.60474308300395
# >CL00000002|90aca794c7e30b8a77e87f13ffc9a5cc;samples=21372;relabund_sum=40.4904515450;size=2629276|100.0

cl_var_seqs = {}
clVarStats = {}

titleRead = False
for line in openfile(cl_vars_file, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()

        title_parts = title.split('|')
        clName = title_parts[0]

        # Example:
        # CL00000001|md5;samples=25791;relabund_sum=77.2086052524;size=5579362|100.0
        md5_and_stats = title_parts[1]
        md5_stats_parts = md5_and_stats.split(';')

        relabund_sum = 0.0
        variant_size = 0

        for item in md5_stats_parts:
            if item.startswith('relabund_sum='):
                relabund_sum = float(item.split('=')[1])
            elif item.startswith('size='):
                variant_size = int(item.split('=')[1])

        if clName not in clVarStats:
            clVarStats[clName] = {
                'relabund_sum': 0.0,
                'variants': 0,
                'variant_size_sum': 0
            }

        clVarStats[clName]['relabund_sum'] += relabund_sum
        clVarStats[clName]['variants'] += 1
        clVarStats[clName]['variant_size_sum'] += variant_size

    else:
        if titleRead:
            titleRead = False
            title_parts = title.split('|')
            clName = title_parts[0]
            seq = line.strip()
            cl_var_seqs[seq] = clName

print('Clustered variants loaded...')

######################
# 2. SAMPLES-STUDIES #
######################

# GB01018568S     Cottin_2025_1acn_Bact
# GB01002033S     Meszarosova_2024_BVF_Bact
# GB01004395S     Netherway_2024_BUB_Bact

sam_stud = {}
for line in openfile(samplesAndStud, 'r'):
    parts = line.strip().split('\t')
    sam_stud[parts[0]] = parts[1]
print("Samples and study names loaded: " + str(len(sam_stud)))

#####################
# 4. REPRESENTATIVE #
#####################

# clName  md5     totSeqSize      totRepSize      HIT     TAXONOMY        SIMILARITY      COVERAGE        EVALUE  BITSCORE        SEQ
# clName  md5     repSamples      repRelAbund_sum repSeqs HIT     TAXONOMY        SIMILARITY      COVERAGE        EVALUE  BITSCORE        SEQ
clRep = {}
with openfile(processed_blast, "r") as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        seq = parts[11]
        clRep[seq] = parts[0]
print('Rep variants loaded...')

#####################
# 3. SAMPLES COUNTS #
#####################

# >GB01020442S|An_2019_1acp_Bact|SRR5920425.11|POS=3|POS=253
# TACAGAGGTCCCAAGCGTTGTTCGGATTTACTGGGCGTAAAGGGCGCGTAGGCGGTTTAGCAAGTTAGAGGTGAAAGGCCCGGGCTTAACCTGGGAACTGCCTTTAAGACTGCTAGGCTTGAGTTCGGAAGAGGATAGCGGAATTCCTAGTGTAGAGGTGAAATTCGTAGATATTAGGAAGAACACCAGTGGCGAAGGCGGCTATCTGGTCCGAAACTGACGCTGAAGCGCGACAGCGTGGGGAGCGAACGGG

clCounts = {}
clStats = {}
clRepStats = {}
titleRead = False
for line in openfile(all_seqs_fasta, 'r'):
    ch = line[0]
    if ch == '>':
        titleRead = True
        sampleName = line[1:].strip().split('|')[0]
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            if cl_var_seqs.has_key(seq):
                clName = cl_var_seqs[seq]
                studName = sam_stud[sampleName]
                if clStats.has_key(clName):
                    clSamp = clStats[clName]['samples']
                    clStud = clStats[clName]['studies']
                    clSamp.add(sampleName)                    
                    clStud.add(studName)
                    clStats[clName]['samples'] = clSamp
                    clStats[clName]['studies'] = clStud
                    # counts
                    clCounts[clName] += 1
                else:
                    clStats[clName] = {}
                    clSamp = set()
                    clStud = set()
                    clSamp.add(sampleName)                    
                    clStud.add(studName)
                    clStats[clName]['samples'] = clSamp
                    clStats[clName]['studies'] = clStud
                    # counts
                    clCounts[clName] = 1
                # representative
                if clRep.has_key(seq):
                    if clRepStats.has_key(clName):
                        clSamp = clRepStats[clName]['samples']
                        clStud = clRepStats[clName]['studies']
                        clSamp.add(sampleName)                    
                        clStud.add(studName)
                        clRepStats[clName]['samples'] = clSamp
                        clRepStats[clName]['studies'] = clStud
                    else:
                        clRepStats[clName] = {}
                        clSamp = set()
                        clStud = set()
                        clSamp.add(sampleName)                    
                        clStud.add(studName)
                        clRepStats[clName]['samples'] = clSamp
                        clRepStats[clName]['studies'] = clStud
            else:
                print("error seq not found "+seq)

print('Sample and study counts were loaded...')

##########################
# 3. REP PROCESSED BLAST #
##########################

# clName  md5     totSeqSize      totRepSize      HIT     TAXONOMY        SIMILARITY      COVERAGE        EVALUE  BITSCORE        SEQ
# CL00000001      6472eb8b1e09f892aca2f23182962903        27490562        6712798 GB-GCA-004799405.1-SSMW01000135.1       d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rhizobiales_505101;f__Xanthobacteraceae;g__Bradyrhizobium_503372;s__Bradyrhizobiumsp000244915   100.000 100.0   1.18e-130       468
#     TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAGGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGG

fp = open(output, 'w')
i = 0

for line in openfile(processed_blast, 'r'):
    line = line.strip()
    parts = line.split("\t")
    clName = parts[0]

    if i == 0:
        fp.write(
            line +
            "\ttotalSeqs"
            "\ttotalSamples"
            "\ttotalStudies"
            "\trepSamples"
            "\trepStudies"
            "\tclusterRelAbund_sum"
            "\tclusterVariants"
            "\tclusterVariantSize_sum"
            "\n"
        )
    else:
        fp.write(
            line + "\t" +
            str(clCounts[clName]) + "\t" +
            str(len(clStats[clName]['samples'])) + "\t" +
            str(len(clStats[clName]['studies'])) + "\t" +
            str(len(clRepStats[clName]['samples'])) + "\t" +
            str(len(clRepStats[clName]['studies'])) + "\t" +
            str(clVarStats[clName]['relabund_sum']) + "\t" +
            str(clVarStats[clName]['variants']) + "\t" +
            str(clVarStats[clName]['variant_size_sum']) +
            "\n"
        )

    i = i + 1

fp.close()

print("DONE!")