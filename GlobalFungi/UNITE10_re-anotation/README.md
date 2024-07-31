### BLAST

- make database

`makeblastdb -in UNITEv10_sh_99_s_all.fasta -dbtype 'nucl' -out UNITEv10_sh_99_s_all`

- run blast

```
for file in *.fas.gz
do  
 echo "gzip -dc ${file} | blastn -query - -db UNITEv10_sh_99_s_all -evalue 1E-5 -outfmt 6 -num_threads 1 -max_target_seqs 10 | gzip > ${file%%.fas.gz}.UNITEv10_sh_99_s_all.txt.gz"
done > blast_command.sh
```


`cat blast_command.sh | parallel`


- get the best hits

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

```
for file in *.txt.gz
do  
  zcat "${file}" | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge | gzip > "${file%%.txt.gz}_best.tab.gz"
done
```

`cat *_best.tab.gz > GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITEv10_sh_99_s_all_BEST_ALL.gz`


### get processed blast file

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/PROCESS_BLAST_RESULTS_SIMPLE.py`

`python2.7 PROCESS_BLAST_RESULTS_SIMPLE.py UNITEv10_sh_99_s_all.fasta GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITEv10_sh_99_s_all_BEST_ALL.gz GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITEv10_sh_99_s_all_PROCESSED.txt`


### GET UNITE 10 All eukaryotes sh dynamic 04042024 complete taxonomy

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/unite10_complete_taxonomy.txt`**


### UPDATE ORIGINAL TABLE

BLASTn results processed - GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz

| QUERY                          | HIT             | SIMILARITY | COVERAGE | EVALUE     | BITSCORE |
|--------------------------------|-----------------|------------|----------|------------|----------|
| 3d9468560d2615a3ee0c1a00c1a143b7 | SH0953900.10FU  | 98.718     | 100.0    | 6.35e-116  | 416      |
| ba029c6222f7dafd1846ff692af82947 | SH0980029.10FU  | 96.622     | 100.0    | 1.84e-64   | 244      |
| 050932ce68b26bfed574a42f9b66cd85 | NO_HIT          | -          | -        | -          | -        |
| 986a02bb45ae7c4fee4d5e78797e3252 | SH0982549.10FU  | 96.266     | 100.0    | 8.54e-110  | 396      |



**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/update_table_processed.py`**

`python2.7 update_table_processed.py GF5_RAW_TABLE_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_VARIANTS_UNITE10_PROCESSED.txt.gz GF5_RAW_TABLE_PROCESSED_UNITE10.txt 98.5 90.0 unite10_complete_taxonomy.txt`


**Output table looks like this...(HEADER IS NOT PART OF THE TABLE)**


| QUERY                          | SAMPLES              | ABUNDANCES      | MARKER | SH       | SEQUENCE |
|--------------------------------|----------------------|-----------------|--------|----------|----------|
| 95e494b4061325c3d39e5f304528191d | 29                   | 1               | ITS1   | 0        | AAAAA... |
| 0fa49bed213860e0db1c33527cfd098f | 30;31;32;33;34;35;36;37;38;39 | 1;1;1;1;1;1;3;1;1;1 | ITS1   | 5        | CCGAG... |
| d2a2561d72610a4dea34c77bd6f3cfc2 | 40                   | 1               | ITS2   | 0        | ACACC... |
| 560dbd59fca8359abdf30da51851ae83 | 41                   | 2               | ITS1   | 0        | CCGAA... |
| 9c783f603a09389a9a2c09aacb3fe250 | 42;43                | 1;1             | ITS2   | 0        | CCACC... |
| df5271d8c7128f9a7e8bbe6e644107da | 44                   | 1               | ITS2   | 6        | AGCCT... |




It generates "UPDATED_TAX_TABLE.txt" containing reduced taxonomy table...


### GET TAXONOMY TABLES

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/get_taxa_table.py`**

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SH.txt 0`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_GENUS.txt 6`

`python2.7 get_taxa_table.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SPECIES.txt 7`


### GET TAXONOMY TABLES FOR SPECIFIC MARKER

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/get_taxa_table_for_marker.py`**

`python2.7 get_taxa_table_for_marker.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SH_ITS1.txt 0 ITS1`

`python2.7 get_taxa_table_for_marker.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt UPDATED_TAX_TABLE.txt GF5_RAW_TABLE_TAB_UNITE10_SH_ITS2.txt 0 ITS2`


### RECALCULATE TABLE FOR distribution modelling

Make abundance table...

**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/OTU_and_SH_tables/TABLE_BY_SAMPLES_AND_SHS.py`**

`python2.7 TABLE_BY_SAMPLES_AND_SHS.py GF5_RAW_TABLE_TAB_UNITE10_SH_ITS1.txt GF5_RAW_TABLE_SAMPLES.txt.gz ELIGIBLE_SAMPLES_ITS1.txt -`


- For each SH, samples where its abundance is > 1 (nonsingletons) AND abundances is > 0.00005 are added to the list where SH is present ("1")
- For each SH, samples where its abundance is = 0 are added to the list where SH is absent ("0")
- For each SH, all other samples (i.e., those where the SH was a local singleton or had very low abundance) are excluded from the sample list ("NA"); it is unclear if their low sequence abundance represents a presence or cross-contamination


Structure of file ELIGIBLE_SAMPLES_ITS1.txt
| PERMANENT_ID | SEQ.COUNTS  |
|--------------|-------------|
| GF01000550S  | 1943        |
| GF01000709S  | 58066       |
| GF01000712S  | 15301       |



**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/calculate_table.py`**


`python2.7 calculate_table.py ELIGIBLE_SAMPLES_ITS1_OUTPUT_TAB.txt ELIGIBLE_SAMPLES_ITS1.txt ELIGIBLE_SAMPLES_ITS1_OUTPUT_DISTR_TAB.txt 0.00005`


### GET STANDALONE VARIANTS FOR TAXA


**`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/UNITE10_re-anotation/GET_VARINATS_FOR_ALL_TAXA_NEW.py`**

`mkdir VARIANTS`

`python2.7 GET_VARINATS_FOR_ALL_TAXA_NEW.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz GF5_UNITE10_TAXONOMY_TABLE.txt.gz VARIANTS/`






python2.7 GET_VARINATS_FOR_ALL_TAXA_FINAL.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz REL3_TAXONOMY.txt VARIANTS/

