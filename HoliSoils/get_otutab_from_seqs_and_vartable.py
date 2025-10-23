__author__ = 'vetrot'

import sys
import gzip
import re

variants = sys.argv[1]  # VARIANTS_variants_holisoils.txt.gz
seqs_sampl = sys.argv[2]  # GB1_samples_holisoils_min10k_max10k.fa.gz
out_tab = sys.argv[3]	# Chao1_Rich_TAB_GB1_10k_samples.txt
taxonomy = sys.argv[4]  # TAXONOMY_CLUSTERS_holisoil.txt

def openfile(filename, mode='r'):
    # For Py3 use 'rt'/'wt' with encoding if needed.
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def get_phylum(tax_string):
    """
    Extract the Phylum name (p__) from a taxonomy string.
    Returns 'unidentified' or None if not found.
    """
    if not tax_string:
        return None

    # Normalize spacing and case
    s = tax_string.strip()

    # Handle completely unidentified cases
    if s.lower() == 'unidentified':
        return 'unidentified'

    # Try to extract using regex
    match = re.search(r'p__([^;]+)', s)
    if match:
        phylum = match.group(1).strip()
        # Handle empty p__ (like p__)
        if phylum == '' or phylum.lower() == 'unidentified':
            return 'unidentified'
        return phylum
    else:
        return 'unidentified'

def chao1_Sobs_per_sample(taxa_abund):
    """
    Compute Chao1 richness estimator and observed species richness for each sample.

    Parameters
    ----------
    taxa_abund : dict
        {sample_id: {taxon_id: count, ...}}  OR  {sample_id: [counts...]}

    Returns
    -------
    dict
        {sample_id: {'S_obs': float, 'Chao1': float}}

    Notes
    -----
    Chao1 = S_obs + F1^2 / (2*F2), if F2 > 0
          = S_obs + F1*(F1-1) / 2,  if F2 == 0 (bias-corrected fallback)
    where:
      S_obs = number of observed taxa with count > 0
      F1    = number of singletons (count == 1)
      F2    = number of doubletons (count == 2)
    """
    out = {}
    for sample, data in taxa_abund.iteritems():
        # Get a flat list of counts for this sample
        if isinstance(data, dict):
            counts = data.values()
        elif isinstance(data, (list, tuple)):
            counts = data
        else:
            # Try to treat as generic iterable
            try:
                counts = list(data)
            except Exception:
                counts = []

        # Filter to positive counts only
        pos = [c for c in counts if c and c > 0]

        S_obs = len(pos)
        F1 = sum(1 for c in pos if c == 1)
        F2 = sum(1 for c in pos if c == 2)

        if S_obs == 0:
            chao1 = 0.0
        else:
            if F2 > 0:
                chao1 = S_obs + (F1 * F1) / (2.0 * F2)
            else:
                # Recommended fallback when no doubletons are observed
                chao1 = S_obs + (F1 * (F1 - 1)) / 2.0

        out[sample] = {
            'S_obs': float(S_obs),
            'Chao1': float(chao1)
        }

    return out



# 8       625     8712a8a6ccf8f2414568cb7cc2f3d92e        TAC...
# 15      281     d3259829740b724d964bf7df4f0b78b2        TAC...

variant_cl = {}
for line in openfile(variants):
	vals = line.rstrip().split('\t')
	cl_id = vals[1]
	seq = vals[3]
	variant_cl[seq] = cl_id

# >GB01015369S|Anthony_2024_CEZ_Bact|SRR27705321.8227|POS=0|POS=253
# TAC...
# >GB01015369S|Anthony_2024_CEZ_Bact|SRR27705321.64886|POS=0|POS=253
# TAC...

titleRead = False
samples_taxa = {}
for line in openfile(seqs_sampl):
    ch = line[0]
    if ch == '>':
        titleRead = True
        title = line[1:].strip()
    else:
        if titleRead:
            titleRead = False
            seq = line.strip()
            sample_name = title.split('|')[0]
            taxa_name = variant_cl[seq]
            if samples_taxa.has_key(sample_name):
            	taxa_abbund = samples_taxa[sample_name]
            	if taxa_abbund.has_key(taxa_name):
            		taxa_abbund[taxa_name] += 1
            	else:
            		taxa_abbund[taxa_name] = 1
            	samples_taxa[sample_name] = taxa_abbund
            else:
            	taxa_abbund = {}
            	taxa_abbund[taxa_name] = 1
            	samples_taxa[sample_name] = taxa_abbund

print("Processed FASTA...")

chao1 = chao1_Sobs_per_sample(samples_taxa)

print("Computed Chao1, Richness...")

# 26      GB00000026.1    -       -       100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rickettsiales; f__Mitochondria; g__; s__     c0d5395792eadbf5f62e8ffb14fa0262
# 1       GB00000001.1    Bradyrhizobium sp000244915      Bradyrhizobium_503372   100.0   100.0   d__Bacteria; p__Pseudomonadota; c__Alphaproteobacteria; o__Rhizobiales_505101; f__Xanthobacteraceae; g__Bradyrhizobium_503372; s__Bradyrhizobium sp000244915    6472eb8b1e09f892aca2f23182962903
# 64      GB00000064.1    -       -       99.605  100.0   d__Bacteria; p__Cyanobacteriota; c__Chloroplast; o__; f__; g__; s__     6427ec8be9d0781c424ff9092dc00dce
# 2       GB00000002.1    Arthrobacter_E_385749 humicola  Arthrobacter_E_385749   100.0   100.0   d__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Actinomycetales; f__Micrococcaceae; g__Arthrobacter_E_385749; s__Arthrobacter_E_385749 humicola   7ff346973a282aa55de296afdb5d74af
# 4147    GB00004147.1    -       -       96.838  100.0   unidentified    00a3dc3cd235b22b6f7497c198fabd5b
allph = set()
cl_phylum = {}
cl_sp_ph = {}
for line in openfile(taxonomy):
    vals = line.rstrip().split('\t')
    cl_id = vals[0]
    sp = vals[2]
    phl = get_phylum(vals[6])
    cl_phylum[cl_id] = phl
    allph.add(phl)
    if not sp == '-':
        if cl_sp_ph.has_key(cl_id):
            phl_sp = cl_sp_ph[cl_id]
            if phl_sp.has_key(sp):
                phl_sp[sp] += 1
            else:
                phl_sp[sp] = 1
            cl_sp_ph[cl_id] = phl_sp
        else:
            phl_sp = {}
            phl_sp[sp] = 1
            cl_sp_ph[cl_id] = phl_sp

phylum_list = []
for phl in allph:
    phylum_list.append(phl)

print("taxonomy loaded... "+str(len(allph)))



phylum_counts = {}
for sample_name in samples_taxa:
    taxa_abbund = samples_taxa[sample_name]
    # get abundances for sample phyla
    phl_abbund = {}
    phl_spp = {}
    for cl_id in taxa_abbund:
        if cl_phylum.has_key(cl_id):
            phl = cl_phylum[cl_id]
            if phl_abbund.has_key(phl):
                phl_abbund[phl] += taxa_abbund[cl_id]
            else:
                phl_abbund[phl] = taxa_abbund[cl_id]
            # sp counts
            if cl_sp_ph.has_key(cl_id):
                phl_sp = cl_sp_ph[cl_id]
                phl_spp[phl] = len(phl_sp)
    # get phylum abbund
    line = ""
    for phl in phylum_list:
        if phl_abbund.has_key(phl):
            line += "\t" + str(phl_abbund[phl])
        else:
            line += "\t0"
    # get phylum sp counts
    for phl in phylum_list:
        if phl_spp.has_key(phl):
            line += "\t" + str(phl_spp[phl])
        else:
            line += "\t0"
    phylum_counts[sample_name] = line

print("phylum counts and sp processed...")


fp = open(out_tab, 'w')
# header
line = "sample\tS_obs\tChao1"
for phl in phylum_list:
    line += "\t" + phl
for phl in phylum_list:
    line += "\t#sp_" + phl
fp.write(line + "\n")
# body
for sample_name in chao1:
	fp.write(sample_name + "\t" + str(chao1[sample_name]['S_obs']) + "\t" + str(chao1[sample_name]['Chao1']) + phylum_counts[sample_name] + '\n')
fp.close()

print("DONE :)")