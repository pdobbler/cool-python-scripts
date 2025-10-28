### PREPARE TABLES FOR FUNGAL SAMPLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/HoliSoils/convert_variants_table_fungi.py`

`python2.7 convert_variants_table_fungi.py GF5_RAW_TABLE_PROCESSED_UNITE10.txt.gz GF5_samples_holisoils.txt`

### SUBSAPLE FASTA BASED ON SAMPLE SIZE RANGE (discard samples under the "discard treshold")

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/get_subsampled_FASTA_FOR_GB.py`

`python get_subsampled_FASTA_FOR_GB.py 10000 GB1_samples_holisoils.fa.gz GB1_samples_holisoils_min10k_max10k.fa GB1_samples_holisoils_discarded_10k.txt 10000`

### GENERATE Chao1 Richness and Phylum breakdown table

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/HoliSoils/get_otutab_from_seqs_and_vartable.py`

`python2.7 get_otutab_from_seqs_and_vartable.py VARIANTS_variants_holisoils.txt.gz GB1_samples_holisoils_min10k_max10k.fa.gz Chao1_Rich_Phylum_TAB_GB1_10k_samples.txt TAXONOMY_CLUSTERS_holisoil.txt`  

### ESTABLISH APP ON EXISTING SERVER

```
cd /home/ubuntu/Docker_holisoils
docker build --network host -t holisoils .
```

edit 
```
proxy:
  port: 8080
  landing-page: /
  heartbeat-rate: 10000
  heartbeat-timeout: 60000
  authentication: none
  hide-navbar: true
specs:
  - id: GlobalFungi
    display-name: GlobalAMFungi
    container-cmd: ["/usr/bin/shiny-server.sh"]
    container-volumes: ["databases:/home/fungal/databases"]
    container-image: fungi_test
  - id: Holisoils
    display-name: Holisoils
    container-image: holisoils
    container-cmd: ["/usr/bin/shiny-server.sh"]
    # separate name of volume:
    container-volumes: ["holisoils_databases:/home/fungal/databases"]
spring:
    servlet:
      multipart:
        max-file-size: 10000MB
        max-request-size: 10000MB

```

### TABLE PREPARATION

- get subsaples papers

```
zcat SAMPLES_PAPERS.txt.gz \
| awk -F'\t' 'NR==FNR { ids[$2]; next } ($1 in ids)' SAMPLES_BASIC_filtered.txt - \
> SAMPLES_PAPERS_finalsamples.txt
```

- get subsaples samplevar data

```
zcat VARIANTS_FUN_samplevar.txt.gz \
| awk -F'\t' 'NR==FNR { ids[$1]; next } ($3 in ids)' SAMPLES_BASIC_filtered.txt - \
> VARIANTS_FUN_samplevar_finalsamples.txt
```

- get varians for selected samples

`awk -F'\t' '!seen[$2]++ {ids[$2]} END {for (i in ids) print i}' VARIANTS_FUN_samplevar_finalsamples.txt | awk 'NR==FNR {ids[$1]; next} $1 in ids' - <(zcat VARIANTS_FUN_variants.txt.gz) > VARIANTS_FUN_variants_finalsamples.txt`

### CREATE mySQL TABLES

```
ALTER TABLE variants_bac
ADD COLUMN marker VARCHAR(4) NOT NULL DEFAULT '16S';
```

```
CREATE TABLE IF NOT EXISTS `variants_fun` (
  `id` int(10) unsigned NOT NULL,
  `cl_id` int(10) unsigned NOT NULL,
  `hash` varchar(32) NOT NULL,
  `sequence` TEXT NOT NULL,
  `marker` VARCHAR(4) NOT NULL
);
```

```
ALTER TABLE variants_bac
  ADD INDEX idx_variants_hash_id_clid (hash, id, cl_id);
ALTER TABLE variants_fun
  ADD INDEX idx_variants_hash_id_clid (hash, id, cl_id);
```
