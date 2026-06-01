### EXTRACT GB region

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/search_for_primary_motive_with_reverse_cut_primary.py`

`python2.7 search_for_primary_motive_with_reverse_cut_primary.py ssu_all_r232.fna GTGYCAGCMGCCGCGGTAA 4`

519F  

```
MOTIVE="GTGYCAGCMGCCGCGGTAA"
MISMATCHES=4

# Run the Python script in parallel for all .fas.gz files
ls *.fas | parallel -j $(nproc) "python2.7 search_for_primary_motive_with_reverse_cut_primary.py {} $MOTIVE $MISMATCHES"
```

806R 

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/find_and_cut_secondary.py`

```
MOTIVE="ATTAGANACCCNNGTA"
MISMATCHES=5

# Run the Python script in parallel for all .PRIMARYCUT.fa.gz files
ls *.PRIMARYCUT.fa.gz | parallel -j $(nproc) "python2.7 find_and_cut_secondary.py {} $MOTIVE $MISMATCHES"
```

### PREPARE TABLES

`zgrep --no-group-separator -wFf ssu_all_final_for_blast_derep_seqs.txt /mnt/DATA1/GLOBAL_BACTERIA/FINAL/DATABASE_TABLES_FINAL/VARIANTS_variants.txt.gz > VARIANTS_variants_GTDB_EXACT.txt`

