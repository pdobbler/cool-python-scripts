### chimera clean

`vsearch --uchime_denovo input.fa.gz --nonchimeras nonchimeras.fa --chimeras chimeras.fa`

### INPUT DATA

/mnt/DATA/projects/avetrot/REL4_RAW_and_CLUSTERING_FINAL/RAW_CLUSTERED/FINAL/REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz

### MAKE VARIANTS WITH SIZE INFO

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz` 

### SEPARATE SINGLE AND MULTIVARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_multi_variants.py`

`python2.7 get_multi_variants.py GF5_ALL_SAMPLES.fa.its1.gz.uniq.gz` 

### CLUSTERING

`vsearch --cluster_size input.uniq.multi --sizein -id 0.97 -uc clusters_97.uc`

### GET CLUSTERS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_vsearch/add_cluster_names_from_vsearch_uc_file.py`

`python2.7 add_cluster_names_from_vsearch_uc_file.py clusters_97.uc REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz.uniq.multi`

### SEEDs IDENTIFICATION

`blastn -query REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz.uniq.multi.seeds -db /mnt/DATA/DATABASES/UNITE10/UNITE_10_SEED2_VERSION_04042024/UNITE10ECOLOGY -out REL4_ITS2_VSEARCH98_unite10.txt -evalue 1E-5 -outfmt 6 -num_threads 128 -max_target_seqs 10`

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr REL4_ITS2_VSEARCH98_unite10.txt | sort -u -k1,1 --merge > REL4_ITS2_VSEARCH98_unite10_best.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`

`python2.7 PROCESS_BLAST_RESULT.py REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz.uniq.multi.seeds REL4_ITS2_VSEARCH98_unite10_best.txt REL4_ITS2_VSEARCH98_unite10_PROCESSED.txt ITS2`

e.g.: Get ectomycorrhizal clusters

`grep 'ectomycorrhizal' REL4_ITS2_VSEARCH98_unite10_PROCESSED.txt > REL4_ITS2_VSEARCH98_unite10_EcM_taxonomy.txt`

`awk -F'|' '{print $1}' REL4_ITS2_VSEARCH97_unite10_EcM_taxonomy.txt > VSEARCH97_unite10_EcM_clusters.txt`

`grep --no-group-separator -A 1 -F -f VSEARCH97_unite10_EcM_clusters.txt REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED.fa.gz.uniq.multi.clustered > REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED_EcM.fa`

### GETTING ABUNDANCE TABLE FROM CLUSTERED VARIANTS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_vsearch/get_table_from_clustered_variants_vsearch.py`

`python2.7 get_table_from_clustered_variants_vsearch.py REL4_ITS2_FUNGAL_AND_NOHIT_CLUSTERED_AND_BINNED_EcM.fa /mnt/DATA/projects/avetrot/REL4_RAW_and_CLUSTERING_FINAL/RELEASE4_RAW/REL4_ITS2_COMPLETE_CLEAN_FUNGAL_AND_NOHIT_FINAL_noReplicated.fa.gz`




