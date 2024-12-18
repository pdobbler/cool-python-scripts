### chimera clean

vsearch --uchime_denovo input.fa.gz --nonchimeras nonchimeras.fa --chimeras chimeras.fa

### INPUT DATA

/mnt/DATA/projects/avetrot/REL4_RAW_and_CLUSTERING_FINAL/RAW_CLUSTERED/FINAL/REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz

### MAKE VARIANTS WITH SIZE INFO

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz` 

### SEPARATE SINGLE AND MULTIVARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_multi_variants.py`

`python2.7 get_multi_variants.py GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz` 

### CLUSTERING

vsearch --cluster_size input.fa --sizein -id 0.97 -uc clusters.uc
