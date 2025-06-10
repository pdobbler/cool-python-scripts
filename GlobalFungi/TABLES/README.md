**Subject:** ITS2 Data and Clustering Information  

As mentioned, there are tables in a reduced format:  
- `GF5_ALL_SAMPLES_ITS2_PERM_OTUTAB_REDUCED_FORMATE.txt.gz`  
- `GF5_ITS2_CLUSTERED_UPARSE97_SAMPLES_in_OTUs_REDUCED_FORMATE.txt.gz`  

These are tab-delimited tables with three columns describing samples in clusters (OTUs) and their abundances, for example:  

`CL114442 GF01016371S;GF01016332S;GF03006580S 3;1;5`


### Clustering Details  
**Permanent Clustering:**  
This was performed based on 98.5% similarity, using sequences present in at least 5 samples ("qualified sequences"). The clusters were then sorted by their relative abundance in samples. After clustering, sequences sufficiently similar (similarity + coverage = 198.5%) were assigned to clusters using BLASTn.  
Details are available here:  
[Permanent Clustering Documentation](https://github.com/pdobbler/cool-python-scripts/tree/main/GlobalFungi/PermanentClusters)  

**UPARSE97 Clustering:**  
Sequence variants appearing at least twice in the dataset were clustered at 97% similarity using the USEARCH 11 UPARSE algorithm. Singletons (variants present only once) were later binned (similarity + coverage = 197%) using BLASTn.  

### Representative Sequences  
Representative sequences for the clusters can be found in these files:  
- `GF5_ALL_SAMPLES_its2_PERMANENT_seeds.fa.gz`  
- `GF5_ITS2_CLUSTERS_UPARSE97_MOST_ABUND.fa.gz`  

### Cluster Annotations  
Cluster annotations based on UNITE10 are available here:  
- `GF5_ALL_SAMPLES_ITS2_CLUSTERED_MOST_ABUND_unite10_PROCESSED.txt.gz`  
- `GF5_ALL_SAMPLES_its2_PERMANENT_seeds_UNITE10_PROCESSED.txt.gz`  

### Additional Script  
I’ve included a script for retrieving a classical tab-delimited OTU table from the reduced format, using files that list samples and clusters:  
`GET_OTUTABLE_BY_SAMPLES_AND_SHS.py`  

**Usage:**  
`python2.7 GET_OTUTABLE_BY_SAMPLES_AND_SHS.py GF5_ALL_SAMPLES_ITS2_PERM_OTUTAB_REDUCED_FORMATE.txt.gz samples.txt clusters.txt OUTPUT_TABLE.txt`

You can replace samples or clusters lists with '-' to retrieve all possible values.

### FILTERING FUNGI

`zcat GF5_ALL_SAMPLES_ITS2_CLUSTERED_MOST_ABUND_unite10_PROCESSED.txt.gz | awk 'BEGIN {FS="\t"; OFS="\t"} NR == 1 || ($10 == "k__Fungi" && ($6 <= 1E-50 || ($4 + $5) >= 180))' > FUNGAL_CLUSTERS/GF5_ALL_SAMPLES_ITS2_CLUSTERED_MOST_ABUND_unite10_PROCESSED_FUNGI.txt`

### GET FUNGAL BREAKDOWN

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/TABLES/get_samples_fungal_sequences.py`

`python2.7 get_samples_fungal_sequences.py GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz REL5_REANOT_UNITE10/GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz GF5_RAW_TABLE_SAMPLES.txt.gz`

### GET TAXA GROUPS BREAKDOWN

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/TABLES/count_taxon_groups_from_table.py`

Example of taxonomy_pairs input file [taxa_pairs_ITS1_ITS2.txt]  
SH0885157.10FU	Fungi_Basidiomycota  
SH0885185.10FU	Fungi_Chytridiomycota  
SH0880802.10FU	Fungi_Ascomycota  
SH0881048.10FU	Fungi_Ascomycota  
...  
SH1958183.10FU	Fungi_Ascomycota  
all_other_fungal_sequences	Fungi_other  
  

`python2.7 count_taxon_groups_from_table.py GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz taxa_pairs_ITS1_ITS2.txt ITS1_ITS2_phyla_counts.txt`











