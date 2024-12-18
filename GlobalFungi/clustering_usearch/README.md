### MAKE VARIANTS WITH SIZE INFO

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py GF5_ALL_SAMPLES.fa.its1.gz` 

### SEPARATE SINGLE AND MULTIVARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_multi_variants.py`

`python2.7 get_multi_variants.py GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz` 

### CLUSTERING

`usearch -cluster_otus GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz.multi -minsize 1 -otus GF5_ALL_SAMPLES.fa.its1_uniq_multi_otus.fa -relabel Otu -uparseout GF5_ALL_SAMPLES_its1_minsize2_uparse.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/PROCESS_UPARSE_v11.0.667_RESULTS.py`

`python2.7 PROCESS_UPARSE_v11.0.667_RESULTS.py GF5_ALL_SAMPLES.fa.its1.gz GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED.fa.gz GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz.multi GF5_ALL_SAMPLES_its1_minsize2_uparse.txt`

Sequences were renamed based on the generated OTUs - total: 1264305616  
Clustered to OTUs     : 1067936446  
Removed chimeras      : 16798002  
Unclustered singletons: 179571168  

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/GET_OTUS_MOST_ABUND_SEQUNCES.py`

`python2.7 GET_OTUS_MOST_ABUND_SEQUNCES.py GF5_ALL_SAMPLES_ITS1_minsize2_CLUSTERED.fa.gz GF5_ALL_SAMPLES_ITS1_CLUSTERED_MOST_ABUND.fa`

```
mkdir SPLIT
makeblastdb -in GF5_ALL_SAMPLES_ITS1_CLUSTERED_MOST_ABUND.fa -dbtype 'nucl' -out SPLIT/ITS1_CLUSTERS
```
