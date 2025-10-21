### PREPARE TABLES FOR FUNGAL SAMPLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/HoliSoils/convert_variants_table_fungi.py`

`python2.7 convert_variants_table_fungi.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz GF5_samples_holisoils.txt`

### SUBSAPLE FASTA BASED ON SAMPLE SIZE RANGE (discard samples under the "discard treshold")

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_subsampled_FASTA_FOR_GB.py`

`python get_subsampled_FASTA_FOR_GB.py 10000 GB1_samples_holisoils.fa.gz GB1_samples_holisoils_min10k_max10k.fa GB1_samples_holisoils_discarded_10k.txt 10000`
