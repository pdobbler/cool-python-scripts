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
| awk -F'\t' 'NR==FNR { ids[$2]; next } ($1 in ids)' _filtered.txt - \
> SAMPLES_PAPERS_finalsamples.txt
```

- get subsaples samplevar data

```
zcat VARIANTS_FUN_samplevar.txt.gz \
| awk -F'\t' 'NR==FNR { ids[$1]; next } ($3 in ids)' _filtered.txt - \
> VARIANTS_FUN_samplevar_finalsamples.txt
```

- get varians for selected samples

`awk -F'\t' '!seen[$2]++ {ids[$2]} END {for (i in ids) print i}' VARIANTS_FUN_samplevar_finalsamples.txt | awk 'NR==FNR {ids[$1]; next} $1 in ids' - <(zcat VARIANTS_FUN_variants.txt.gz) > VARIANTS_FUN_variants_finalsamples.txt`

- get taxonomy for selected variants

```
zcat GF5_UNITE10_TAXONOMY_TABLE_processed.txt.gz \
| awk -F'\t' 'NR==FNR { ids[$2]; next } ($1 in ids)' FUN_VARIANTS_variants_finalsamples.txt - \
> GF5_UNITE10_TAXONOMY_TABLE_finalsamples.txt
```

### PREPARE DATABASE

`/home/ubuntu/mysql-data/HOLISOILS`
`scp -i /mnt/DATA1/KEYS/FMT_2.pem FUN_* ubuntu@XXX.XXX.XXX.XXX:/home/ubuntu/mysql-data/HOLISOILS`

`docker exec -it mariadb_ok mariadb -u root -p`
[root]

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

```
UPDATE clusters_tax_fun
SET hash = '-';
```

### FILL DATABASE

```
CREATE TABLE `info` (
  `id` int(10) UNSIGNED NOT NULL,
  `name` text NOT NULL,
  `version` text NOT NULL,
  `release` text NOT NULL,
  `annotation_version` text NOT NULL,
  `variants_count` int(11) NOT NULL,
  `raw_count` int(11) NOT NULL,
  `info` text NOT NULL,
  `citation` text NOT NULL,
  `date` varchar(32) NOT NULL
);
```

ALTER TABLE info CONVERT TO CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
SET NAMES 'utf8mb4' COLLATE 'utf8mb4_unicode_ci';

```
INSERT INTO `info` (`id`, `name`, `version`, `release`, `annotation_version`, `variants_count`, `raw_count`, `info`, `citation`, `date`) VALUES
(1, 'Forestmicrobes', 'v1.0', '1', '(31.10.2025)', 51603682, 1116074647, 'First release', 'Větrovský T., Martinovič T., Baldrian P.: Forestmicrobes database.', '31.10.2025');
```

```
CREATE TABLE `traffic` (
  `id` int NOT NULL PRIMARY KEY,
  `session` int(11) NOT NULL,
  `category` varchar(32) DEFAULT NULL,
  `value` varchar(64) DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT current_timestamp()
);
```

```
CREATE TABLE IF NOT EXISTS `maillist` (
  `id` int unsigned NOT NULL auto_increment PRIMARY KEY,
  `name` TEXT NOT NULL,
  `email` TEXT NOT NULL,
  `date` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);
```

```
CREATE TABLE IF NOT EXISTS `messages` (
  `id` int unsigned NOT NULL auto_increment PRIMARY KEY,
  `email` TEXT NOT NULL,
  `subject` TEXT NOT NULL,
  `message` TEXT NOT NULL,
  `processed` boolean not null default 0,
  `date` TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);
```
  
- SAMPLES BASIC
  
```
CREATE TABLE IF NOT EXISTS `samples_basic_fun` (
  `id` int NOT NULL PRIMARY KEY,
  `paper` int NOT NULL,
  `permanent_id` VARCHAR(16) NOT NULL,
  `sample_type` VARCHAR(32) NOT NULL,
  `latitude` float NOT NULL,
  `longitude` float NOT NULL,
  `continent` VARCHAR(14) NOT NULL,
  `year_of_sampling_from` int,
  `year_of_sampling_to` int,
  `Biome` VARCHAR(32) NOT NULL,
  `primers` VARCHAR(128) NOT NULL,
  `MAT` FLOAT,
  `MAP` FLOAT,
  `pH` FLOAT,
  `SOC` FLOAT,
  `ITS_total` int NOT NULL,
  `manipulated` TINYINT(1) NOT NULL
);
```
```
LOAD DATA INFILE '/var/lib/mysql/HOLISOILS/FUN_SAMPLES_BASIC_filtered.txt'
INTO TABLE `samples_basic_fun`
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(id, paper, permanent_id, sample_type, latitude, longitude, continent, @year_of_sampling_from, @year_of_sampling_to, Biome, primers, @MAT, @MAP, @pH, @SOC, ITS_total, manipulated)
SET
  year_of_sampling_from = NULLIF(@year_of_sampling_from, 'NULL'),
  year_of_sampling_to = NULLIF(@year_of_sampling_to, 'NULL'),
  MAT = NULLIF(@MAT, 'NULL'),
  MAP = NULLIF(@MAP, 'NULL'),
  pH = NULLIF(@pH, 'NULL'),
  SOC = NULLIF(@SOC, 'NULL');
```

`ALTER TABLE samples_basic_fun ADD INDEX(id);`

`ALTER TABLE samples_basic_fun CHANGE ITS_total seqs_total INT;`

```
CREATE TABLE IF NOT EXISTS `samples_basic_bac` (
  `id` int NOT NULL PRIMARY KEY,
  `paper` int NOT NULL,
  `permanent_id` VARCHAR(16) NOT NULL,
  `sample_type` VARCHAR(32) NOT NULL,
  `latitude` float NOT NULL,
  `longitude` float NOT NULL,
  `continent` VARCHAR(14) NOT NULL,
  `year_of_sampling_from` int,
  `year_of_sampling_to` int,
  `Biome` VARCHAR(32) NOT NULL,
  `primers` VARCHAR(128) NOT NULL,
  `MAT` FLOAT,
  `MAP` FLOAT,
  `pH` FLOAT,
  `SOC` FLOAT,
  `ITS_total` int NOT NULL,
  `manipulated` TINYINT(1) NOT NULL
);
```
```
LOAD DATA INFILE '/var/lib/mysql/HOLISOILS/BAC_SAMPLES_BASIC_final.txt'
INTO TABLE `samples_basic_bac`
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(id, paper, permanent_id, sample_type, latitude, longitude, continent, @year_of_sampling_from, @year_of_sampling_to, Biome, primers, @MAT, @MAP, @pH, @SOC, ITS_total, manipulated)
SET
  year_of_sampling_from = NULLIF(@year_of_sampling_from, 'NULL'),
  year_of_sampling_to = NULLIF(@year_of_sampling_to, 'NULL'),
  MAT = NULLIF(@MAT, 'NULL'),
  MAP = NULLIF(@MAP, 'NULL'),
  pH = NULLIF(@pH, 'NULL'),
  SOC = NULLIF(@SOC, 'NULL');
```

`ALTER TABLE samples_basic_bac ADD INDEX(id);`

`ALTER TABLE samples_basic_bac CHANGE ITS_total seqs_total INT;`

- SAMPLES ADVANCED

```
CREATE TABLE IF NOT EXISTS `samples_advanced_fun` (
  `id` int NOT NULL PRIMARY KEY,
  `sample_name` VARCHAR(128) NOT NULL,
  `sample_description` TEXT NOT NULL,
  `sequencing_platform` VARCHAR(16) NOT NULL,
  `target_gene` VARCHAR(7) NOT NULL,
  `primers_sequence` VARCHAR(256) NOT NULL,
  `sample_seqid` VARCHAR(256) NOT NULL,
  `sample_barcode` VARCHAR(128) NOT NULL,
  `elevation` INT,
  `MAT_study` FLOAT,
  `MAP_study` FLOAT,
  `Biome_detail` VARCHAR(64) NOT NULL,
  `country` VARCHAR(64) NOT NULL,
  `month_of_sampling` VARCHAR(32) NOT NULL,
  `day_of_sampling` VARCHAR(16) NOT NULL,
  `plants_dominant` TEXT NOT NULL,
  `plants_all` TEXT NOT NULL,
  `area_sampled` FLOAT,
  `number_of_subsamples_from` INT,
  `number_of_subsamples_to` INT,
  `sampling_info` TEXT NOT NULL,
  `sample_depth_from` FLOAT,
  `sample_depth_to` FLOAT,
  `extraction_DNA_mass_from` FLOAT,
  `extraction_DNA_mass_to` FLOAT,
  `extraction_DNA_size` VARCHAR(256) NOT NULL,
  `extraction_DNA_method` VARCHAR(512) NOT NULL,
  `total_C_content` FLOAT,
  `total_N_content` FLOAT,
  `organic_matter_content` FLOAT,
  `pH_study` FLOAT,
  `pH_method` VARCHAR(12) NOT NULL,
  `total_Ca` FLOAT,
  `total_P` FLOAT,
  `total_K` FLOAT,
  `sample_info` TEXT NOT NULL,
  `location` VARCHAR(256) NOT NULL,
  `area_GPS_from` FLOAT,
  `area_GPS_to` FLOAT,
  `ITS1_extracted` INT NOT NULL,
  `ITS2_extracted` INT NOT NULL
);
```

```
LOAD DATA INFILE '/var/lib/mysql/HOLISOILS/FUN_SAMPLES_ADVANCED_filtered.txt'
INTO TABLE `samples_advanced_fun`
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(
  `id`, `sample_name`, `sample_description`, `sequencing_platform`, `target_gene`, 
  `primers_sequence`, `sample_seqid`, `sample_barcode`, @elevation, @MAT_study, 
  @MAP_study, `Biome_detail`, `country`, `month_of_sampling`, `day_of_sampling`, 
  `plants_dominant`, `plants_all`, @area_sampled, @number_of_subsamples_from, 
  @number_of_subsamples_to, sampling_info, @sample_depth_from, @sample_depth_to, 
  @extraction_DNA_mass_from, @extraction_DNA_mass_to, `extraction_DNA_size`, 
  `extraction_DNA_method`, @total_C_content, @total_N_content, @organic_matter_content, 
  @pH_study, `pH_method`, @total_Ca, @total_P, @total_K, sample_info, `location`, 
  @area_GPS_from, @area_GPS_to, `ITS1_extracted`, `ITS2_extracted`
)
SET
`elevation` = NULLIF(@elevation, 'NULL'),
`MAT_study` = NULLIF(@MAT_study, 'NULL'),
`MAP_study` = NULLIF(@MAP_study, 'NULL'),
`area_sampled` = NULLIF(@area_sampled, 'NULL'),
`number_of_subsamples_from` = NULLIF(@number_of_subsamples_from, 'NULL'),
`number_of_subsamples_to` = NULLIF(@number_of_subsamples_to, 'NULL'),
`sample_depth_from` = NULLIF(@sample_depth_from, 'NULL'),
`sample_depth_to` = NULLIF(@sample_depth_to, 'NULL'),
`extraction_DNA_mass_from` = NULLIF(@extraction_DNA_mass_from, 'NULL'),
`extraction_DNA_mass_to` = NULLIF(@extraction_DNA_mass_to, 'NULL'),
`total_C_content` = NULLIF(@total_C_content, 'NULL'),
`total_N_content` = NULLIF(@total_N_content, 'NULL'),
`organic_matter_content` = NULLIF(@organic_matter_content, 'NULL'),
`pH_study` = NULLIF(@pH_study, 'NULL'),
`total_Ca` = NULLIF(@total_Ca, 'NULL'),
`total_P` = NULLIF(@total_P, 'NULL'),
`total_K` = NULLIF(@total_K, 'NULL'),
`area_GPS_from` = NULLIF(@area_GPS_from, 'NULL'),
`area_GPS_to` = NULLIF(@area_GPS_to, 'NULL');
```

`ALTER TABLE samples_advanced_fun ADD INDEX(id);`

```
CREATE TABLE IF NOT EXISTS `samples_advanced_bac` (
  `id` int NOT NULL PRIMARY KEY,
  `sample_name` VARCHAR(128) NOT NULL,
  `sample_description` TEXT NOT NULL,
  `sequencing_platform` VARCHAR(16) NOT NULL,
  `target_gene` VARCHAR(7) NOT NULL,
  `primers_sequence` VARCHAR(256) NOT NULL,
  `sample_seqid` VARCHAR(256) NOT NULL,
  `sample_barcode` VARCHAR(128) NOT NULL,
  `elevation` INT,
  `MAT_study` FLOAT,
  `MAP_study` FLOAT,
  `Biome_detail` VARCHAR(64) NOT NULL,
  `country` VARCHAR(64) NOT NULL,
  `month_of_sampling` VARCHAR(32) NOT NULL,
  `day_of_sampling` VARCHAR(16) NOT NULL,
  `plants_dominant` TEXT NOT NULL,
  `plants_all` TEXT NOT NULL,
  `area_sampled` FLOAT,
  `number_of_subsamples_from` INT,
  `number_of_subsamples_to` INT,
  `sampling_info` TEXT NOT NULL,
  `sample_depth_from` FLOAT,
  `sample_depth_to` FLOAT,
  `extraction_DNA_mass_from` FLOAT,
  `extraction_DNA_mass_to` FLOAT,
  `extraction_DNA_size` VARCHAR(256) NOT NULL,
  `extraction_DNA_method` VARCHAR(512) NOT NULL,
  `total_C_content` FLOAT,
  `total_N_content` FLOAT,
  `organic_matter_content` FLOAT,
  `pH_study` FLOAT,
  `pH_method` VARCHAR(12) NOT NULL,
  `total_Ca` FLOAT,
  `total_P` FLOAT,
  `total_K` FLOAT,
  `sample_info` TEXT NOT NULL,
  `location` VARCHAR(256) NOT NULL,
  `area_GPS_from` FLOAT,
  `area_GPS_to` FLOAT,
  `ITS1_extracted` INT NOT NULL,
  `ITS2_extracted` INT NOT NULL
);
```

```
LOAD DATA INFILE '/var/lib/mysql/HOLISOILS/BAC_SAMPLES_ADVANCED_final.txt'
INTO TABLE `samples_advanced_bac`
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
(
  `id`, `sample_name`, `sample_description`, `sequencing_platform`, `target_gene`, 
  `primers_sequence`, `sample_seqid`, `sample_barcode`, @elevation, @MAT_study, 
  @MAP_study, `Biome_detail`, `country`, `month_of_sampling`, `day_of_sampling`, 
  `plants_dominant`, `plants_all`, @area_sampled, @number_of_subsamples_from, 
  @number_of_subsamples_to, sampling_info, @sample_depth_from, @sample_depth_to, 
  @extraction_DNA_mass_from, @extraction_DNA_mass_to, `extraction_DNA_size`, 
  `extraction_DNA_method`, @total_C_content, @total_N_content, @organic_matter_content, 
  @pH_study, `pH_method`, @total_Ca, @total_P, @total_K, sample_info, `location`, 
  @area_GPS_from, @area_GPS_to, `ITS1_extracted`, `ITS2_extracted`
)
SET
`elevation` = NULLIF(@elevation, 'NULL'),
`MAT_study` = NULLIF(@MAT_study, 'NULL'),
`MAP_study` = NULLIF(@MAP_study, 'NULL'),
`area_sampled` = NULLIF(@area_sampled, 'NULL'),
`number_of_subsamples_from` = NULLIF(@number_of_subsamples_from, 'NULL'),
`number_of_subsamples_to` = NULLIF(@number_of_subsamples_to, 'NULL'),
`sample_depth_from` = NULLIF(@sample_depth_from, 'NULL'),
`sample_depth_to` = NULLIF(@sample_depth_to, 'NULL'),
`extraction_DNA_mass_from` = NULLIF(@extraction_DNA_mass_from, 'NULL'),
`extraction_DNA_mass_to` = NULLIF(@extraction_DNA_mass_to, 'NULL'),
`total_C_content` = NULLIF(@total_C_content, 'NULL'),
`total_N_content` = NULLIF(@total_N_content, 'NULL'),
`organic_matter_content` = NULLIF(@organic_matter_content, 'NULL'),
`pH_study` = NULLIF(@pH_study, 'NULL'),
`total_Ca` = NULLIF(@total_Ca, 'NULL'),
`total_P` = NULLIF(@total_P, 'NULL'),
`total_K` = NULLIF(@total_K, 'NULL'),
`area_GPS_from` = NULLIF(@area_GPS_from, 'NULL'),
`area_GPS_to` = NULLIF(@area_GPS_to, 'NULL');
```

`ALTER TABLE samples_advanced_bac ADD INDEX(id);`

- PAPERS

```
CREATE TABLE IF NOT EXISTS `samples_papers_fun` (
  `id` int NOT NULL PRIMARY KEY,
  `add_date` VARCHAR(10) NOT NULL,
  `year` int NOT NULL,
  `title` TEXT NOT NULL,
  `authors` TEXT NOT NULL,
  `journal` VARCHAR(128) NOT NULL,
  `doi` VARCHAR(64) NOT NULL,
  `contact` TEXT NOT NULL,
  `manipulated` TINYINT(1) NOT NULL
);  
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/HOLISOILS/FUN_SAMPLES_PAPERS_finalsamples.txt' INTO TABLE samples_papers_fun FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

`ALTER TABLE samples_papers_fun ADD INDEX(id);`

```
CREATE TABLE IF NOT EXISTS `samples_papers_bac` (
  `id` int NOT NULL PRIMARY KEY,
  `add_date` VARCHAR(10) NOT NULL,
  `year` int NOT NULL,
  `title` TEXT NOT NULL,
  `authors` TEXT NOT NULL,
  `journal` VARCHAR(128) NOT NULL,
  `doi` VARCHAR(64) NOT NULL,
  `contact` TEXT NOT NULL,
  `manipulated` TINYINT(1) NOT NULL
);  
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/HOLISOILS/BAC_SAMPLES_PAPERS_final.txt' INTO TABLE samples_papers_bac FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

`ALTER TABLE samples_papers_bac ADD INDEX(id);`

- TAXONOMY

```
CREATE TABLE IF NOT EXISTS `clusters_tax_fun` (
  `id` int NOT NULL PRIMARY KEY,
  `cluster` VARCHAR(12) NOT NULL,
  `Species` varchar(64) NOT NULL,
  `Genus` varchar(32) NOT NULL,
  `sim` FLOAT,
  `cov` FLOAT,
  `full_tax` TEXT NOT NULL,
  `hash` varchar(32) NOT NULL
);
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/HOLISOILS/FUN_UNITE10_TAXONOMY_TABLE_finalsamples.txt' INTO TABLE clusters_tax_fun FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

```
ALTER TABLE clusters_tax_fun ADD INDEX(id);
ALTER TABLE clusters_tax_fun ADD INDEX(cluster);
ALTER TABLE clusters_tax_fun ADD INDEX(Species);
ALTER TABLE clusters_tax_fun ADD INDEX(Genus);
CREATE INDEX idx_clusters_tax_species_id ON clusters_tax_fun (Species, id);
CREATE INDEX idx_clusters_tax_genus_id ON clusters_tax_fun (Genus, id);
```

```
CREATE TABLE IF NOT EXISTS `clusters_tax_bac` (
  `id` int NOT NULL PRIMARY KEY,
  `cluster` VARCHAR(12) NOT NULL,
  `Species` varchar(64) NOT NULL,
  `Genus` varchar(32) NOT NULL,
  `sim` FLOAT,
  `cov` FLOAT,
  `full_tax` TEXT NOT NULL,
  `hash` varchar(32) NOT NULL
);
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/HOLISOILS/BAC_TAXONOMY_CLUSTERS_finalsamples.txt' INTO TABLE clusters_tax_bac FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

```
ALTER TABLE clusters_tax_bac ADD INDEX(id);
ALTER TABLE clusters_tax_bac ADD INDEX(cluster);
ALTER TABLE clusters_tax_bac ADD INDEX(Species);
ALTER TABLE clusters_tax_bac ADD INDEX(Genus);
CREATE INDEX idx_clusters_tax_species_id ON clusters_tax_bac (Species, id);
CREATE INDEX idx_clusters_tax_genus_id ON clusters_tax_bac (Genus, id);
```

- VARIANTS

```
CREATE TABLE IF NOT EXISTS `variants` (
  `id` int(10) unsigned NOT NULL,
  `cl_id` int(10) unsigned NOT NULL,
  `hash` varchar(32) NOT NULL,
  `sequence` TEXT NOT NULL
);
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/GB1_TABLES_RAW/VARIANTS_variants.txt' INTO TABLE variants FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

```
ALTER TABLE variants
  ADD INDEX idx_variants_hash_id_clid (hash, id, cl_id);
```

```
CREATE TABLE IF NOT EXISTS `samplevar` (
  `id` bigint(20) unsigned NOT NULL,
  `variant` int(10) unsigned NOT NULL,
  `sample` int(10) unsigned NOT NULL,
  `abundance` int(10) unsigned NOT NULL,
  `cl_id` int(10) unsigned NOT NULL
);
```

`LOAD DATA LOCAL INFILE '/var/lib/mysql/GB1_TABLES_RAW/VARIANTS_samplevar.txt' INTO TABLE samplevar FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

- taxa search
```
alter table samplevar add index idx_samplevar_clid_sample_abundance (cl_id, sample, abundance);
ALTER TABLE samplevar ADD INDEX idx_samplevar_variant_sample_abundance (variant, sample, abundance);
```
-geosearch
`CREATE INDEX idx_samplevar_sample_clid ON samplevar (sample, cl_id);`

Update stats after creating indexes
`ANALYZE TABLE variants, samplevar, clusters_tax;`
