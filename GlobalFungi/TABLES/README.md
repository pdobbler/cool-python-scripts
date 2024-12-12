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
`python2.7 GET_OTUTABLE_BY_SAMPLES_AND_SHS.py GF5_ALL_SAMPLES_ITS2_PERM_OTUTAB_REDUCED_FORMATE.txt.gz samples.txt clusters.txt`

You can replace samples or clusters lists with '-' to retrieve all possible values.






