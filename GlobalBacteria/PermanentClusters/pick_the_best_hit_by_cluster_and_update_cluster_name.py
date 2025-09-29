__author__ = 'vetrot'

import sys
import os
import hashlib
import gzip
import re

blast_10results = sys.argv[1] # GB_VOL1_20250526_CLEANED_ranked_single3part117_SEEDS97.txt
blast_input_fasta = sys.argv[2] # GB_VOL1_20250526_CLEANED_ranked_single3part117.fas
cluster_names = sys.argv[3]   # /mnt/DATA1/GLOBAL_BACTERIA/FINAL/SEEDS_97.0_WORKING_NAMES.fa.info
treshold = sys.argv[4]  # 197.0 
out_type = sys.argv[5]  # best/top - best = best hit by score, e-val etc.; top - top cluster (smallest) based on treshold 
out_fasta = sys.argv[6] 

if out_type == "best":
    print("Output type: " + out_type + " = best hit by score, e-val etc.")
elif out_type == "top":
    print("Output type: " + out_type + " = top - top cluster (smallest) based on treshold.")
else:
    print("Wrong output type! Quitting...")
    sys.exit(1)

THRESH = float(treshold)

#############################################
# HELPERS
#############################################
def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def hit_numeric_key(hit):
    """
    Extract the smallest integer found in the hit name (e.g. 'CL000624.1' -> 624).
    Falls back to a large number if none found.
    """
    import re
    m = re.search(r'(\d+)(?=\.)', hit)
    if m:
        return int(m.group(1))
    m = re.search(r'(\d+)', hit)
    return int(m.group(1)) if m else 10**12

def parse_float(x):
    return float(x.strip())

def getBestResult(res):
    """
    Pick the typical 'best BLAST hit':
    sort by: bitscore DESC, evalue ASC, pident DESC,
             then (as a stable tie-breaker) by the numeric part of cluster_name ASC.
    """
    if not res:
        return None
    return sorted(
        res,
        key=lambda r: (-r[4], r[3], -r[1], hit_numeric_key(r[0]))
    )[0]


def getTopResult(res):
    """
    From hits passing the threshold (pident + coverage >= THRESH),
    return the best hit (smallest numeric cluster, then bitscore, etc.).
    If ALL hits pass, append an extra 'all_pass' marker to the returned list.
    """
    if not res:
        return None

    ok = [r for r in res if (r[1] + r[2]) >= THRESH]
    if not ok:
        return None

    best = min(
        ok,
        key=lambda r: (hit_numeric_key(r[0]), -r[4], r[3], -r[1])
    )

    # If all hits passed - add "all_pass" marker
    if len(ok) == len(res):
        return best + ["all_pass"]

    return best


#############################################
# LOAD CLUSTER NAMES
#############################################

# GB00000001.1    p__Pseudomonadota|CL000001|6472eb8b1e09f892aca2f23182962903|V_4505598|S_15048|P_97|r_65.6430231571|SEED
# GB00000002.1    p__Actinomycetota|CL000001|7ff346973a282aa55de296afdb5d74af|V_3348910|S_11165|P_91|r_46.9649270943|SEED

clusters = {}
for line in openfile(cluster_names):
    vals = line.rstrip().split('\t')
    clusters[vals[1]] = vals[0]

print("Cluster names loaded...")

# >64999043bd25cd9ad4ef90cf846c5401|V_1|S_1|P_1|r_3.77318622938e-06
# TACGTAGGGGGCTGGCGTTGTCCGGATTTATTGGGCGTACAGCGCGTGTAGGCGGCCGGCTAGGTCTGGTGTGAAAACTCGAGTCTCAACCTCGAGATTTCGCCGGAAACCAGTCGGCTAGAGTCCGGAAGAGGAGAGTGGAATTCCTGGTGTAGCGGTGAACTGCGCAGATCTCAGGAAGAACACCTATGGCGACAGCAGCGCTCTGGGACGGTACTGACGCTGAGCCGCGAAAGCGTGGGGAGCGAACAGG

# load fasta seqs
seqs = {}
titleRead = False
for line in openfile(blast_input_fasta):
    ch = line[0]
    if ch == '>':
        titleRead=True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead=False
            seq = line.strip()
            seqs[title] = seq

print("FASTA loaded...")

# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL000840|ba8f6b576295cab3f871f34ebc59f046|V_2393|S_290|P_34|r_0.0200984047252|SEED   98.419  253     4       0       1       253     1       253     8.68e-124       446
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL338996|2a074b25c9484c0ee56b5e173ea12742|V_2|S_1|P_1|r_4.80136742944e-06|SEED       98.024  253     5       0       1       253     1       253     4.04e-122       440
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL000624|4bf501bbd3ff596ca54f8c949cbc5be6|V_71|S_16|P_6|r_0.035921944756|SEED        97.233  253     7       0       1       253     1       253     8.75e-119       429
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL232355|e72e3ccad093a50716e6d52b332ea45c|V_2|S_1|P_1|r_9.57001904434e-06|SEED       97.222  252     7       0       2       253     2       253     3.15e-118       427
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL343296|fe993bda1c7e9cf78bf62b0e741064fd|V_2|S_1|P_1|r_4.67241685434e-06|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL268416|5b280659f3f7832df89aa82a21eca894|V_2|S_1|P_1|r_7.23541882222e-06|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL229802|82cce0b0d183f5c953ccaf1cd5aa6c81|V_2|S_1|P_1|r_9.70478057491e-06|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL076388|7d2d84af395e3377d8677a7fcab85bcf|V_2|S_1|P_1|r_4.44237133782e-05|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL049524|9377e19a58afc9923e27032c0b2bd5be|V_2|S_2|P_1|r_6.22656239436e-05|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# b4cdc27fe15c034c3cffc2ce9d5ea865|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL027183|9266f24adbeb8bdfc02de51bb80118fe|V_3|S_3|P_1|r_0.000100075504544|SEED       96.838  253     8       0       1       253     1       253     4.07e-117       424
# 8987bb602afafaef5a78c0b6c7f03bc2|V_1|S_1|P_1|r_4.20654117152e-06        p__Acidobacteriota|CL325103|6ce6d7a250dc6556dee751c1a447faca|V_2|S_1|P_1|r_5.35180061332e-06|SEED       99.213  254     2       0       1       254     30      283     1.12e-127       459


# load blast results...
blast_info = {}
for line in openfile(blast_10results):
    values = line.strip().split("\t")
    # storing - cluster_name | pident | coverage |  evalue | bitscore
    slen = len(seqs[values[0]])
    cov = 100.0 * ((int(values[7]) - int(values[6]) + 1) / float(slen))
    val = [clusters[values[1]] ,float(values[2]), float(cov), float(values[10]), float(values[11])]
    if blast_info.has_key(values[0]):
        res = blast_info[values[0]]
        res.append(val)
        blast_info[values[0]] = res
    else:
        res = []
        res.append(val)
        blast_info[values[0]] = res

print("BLAST loaded...")

fp = open(out_fasta + ".all_pass", 'w')
of = gzip.open(out_fasta + ".gz", 'wt')
for seq_name in seqs:
    if blast_info.has_key(seq_name):
        res = blast_info[seq_name]
        result = []
        if out_type == "best":
            result = getBestResult(res)
            of.write('>'+ result[0] +'|' + seq_name + '|BEST_sim' + str(result[1]) + '_cov' + str(result[2]) + '_eval' + str(result[3]) + '_score' + str(result[4]) + '\n')
            of.write(seqs[seq_name] + '\n')
        else:
            result = getTopResult(res)
            if result == None:
                of.write('>NO_HIT|' + seq_name + '|0.0' + '\n')
                of.write(seqs[seq_name] + '\n')
            else:
                if len(result) == 6:
                    of.write('>'+ result[0] +'|' + seq_name + '|' + str(result[1]) + '_all_pass\n')
                    of.write(seqs[seq_name] + '\n')
                    # save seqs to re-blast
                    fp.write('>' + seq_name + '\n')
                    fp.write(seqs[seq_name] + '\n')
                else:
                    of.write('>'+ result[0] +'|' + seq_name + '|' + str(result[1]) + '_ok\n')
                    of.write(seqs[seq_name] + '\n')
    else:
        of.write('>NO_HIT|' + seq_name + '|0.0' + '\n')
        of.write(seqs[seq_name] + '\n')

of.close()
fp.close()

print("DONE :)")

