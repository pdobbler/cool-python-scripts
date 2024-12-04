### MAKE VARIANTS WITH SIZE INFO

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py GF5_ALL_SAMPLES.fa.its1.gz` 

### SEPARATE SINGLE AND MULTIVARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_multi_variants.py`

`python2.7 get_multi_variants.py GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz` 

### CLUSTERING

usearch -cluster_otus GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz.multi -minsize 1 -otus GF5_ALL_SAMPLES.fa.its1_uniq_multi_otus.fa -relabel Otu -uparseout GF5_ALL_SAMPLES_its1_minsize2_uparse.txt

