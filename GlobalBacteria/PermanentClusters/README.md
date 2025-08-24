### GREENGENES TAXONOMY BREAKDOWN

Blasting all variants against Greengenes2 20250505 -> process the blast best hit

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy.zip`  

`unzip taxonomy.zip`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/taxonomy_breakdown.py`  

`python2.7 taxonomy_breakdown.py GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_20250526_CLEANED_ALL_Phylum_breakdown.txt 188.0`

### RANK AND SORT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/generate_SV_ranks_relative_abund.py`

`python2.7 generate_SV_ranks_relative_abund.py GB_VOL1_20250526_CLEANED.fa.gz false`

373 267 631 variants (GB_VOL1_20250526_CLEANED.fa.gz.all)

### FILTER NON-SINGLETONS

`awk '/^>/ {keep = ($0 !~ /\|V_1\|S_1\|/)} keep' GB_VOL1_20250526_CLEANED.fa.gz.all > GB_VOL1_20250526_CLEANED_ranked_multi.fa`

51 603 682 variants (GB_VOL1_20250526_CLEANED_ranked_multi.fa)

`awk '/^>/ {keep = ($0 ~ /\|V_1\|S_1\|/)} keep' GB_VOL1_20250526_CLEANED.fa.gz.all > GB_VOL1_20250526_CLEANED_ranked_single.fa`


### SPLIT FASTA BY TAXONOMY

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/split_by_taxonomy.py`  

`python2.7 split_by_taxonomy.py GB_VOL1_20250526_CLEANED_ranked_multi.fa GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt taxonomy.tsv 1 GB_VOL1_PHYLUM 188.0`

### STATS

```
for file in *.fas; do
    count=$(grep -c '^>' "$file")
    echo -e "$file\t$count" >> counts.txt
done
```

### GROUP SEQUENCES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/group_qualified_sequences_fastest_optimized.py`

```
export TMPDIR=/mnt/DATA1/tmp
parallel --tmpdir /mnt/DATA1/tmp -j $(nproc) python group_qualified_sequences_fastest_optimized.py {} 97.0 false ::: *.fas > output.txt
```

### TRY TO CONTINUE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/split_by_last_seed.py`

`python2.7 split_by_last_seed.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas GB_VOL1_PHYLUM_p__Pseudomonadota.fas.97.0.seed_seqs`

```
mkdir Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.done Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.undone Pseudomonadota_CONTINUE
mv GB_VOL1_PHYLUM_p__Pseudomonadota.fas.input_seeds_seqs Pseudomonadota_CONTINUE
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/group_qualified_sequences_fastest_optimized_seeds.py`

`python group_qualified_sequences_fastest_optimized_seeds.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas.done 97.0 false`

`python group_qualified_sequences_fastest_optimized_seeds.py GB_VOL1_PHYLUM_p__Pseudomonadota.fas.undone 97.0 false GB_VOL1_PHYLUM_p__Pseudomonadota.fas.input_seeds_seqs`


### MULTIPLE - CONTINUE

```
for file in *.fas
do  
 python2.7 split_by_last_seed.py ${file} ${file}.97.0.seed_seqs
done
```

```
for file in *.done
do
  echo "python group_qualified_sequences_fastest_optimized_seeds.py ${file} 97.0 false"
  echo "python group_qualified_sequences_fastest_optimized_seeds.py ${file%%.done}.undone 97.0 false ${file%%.done}.input_seeds_seqs"
done > continue.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat continue.sh | parallel --tmpdir /mnt/DATA1/tmp > output.txt
```

### SEEDS WORKING NAMES

`grep '>' SEEDS_97.0_WORKING_NAMES.fa | wc -l`  
3455524

example:  
>NO_HIT|CL00001|6ccb96db9d81f6f110fe7cb5be5bdf5b|V_15083|S_241|P_11|r_0.143077314152|SEED  
TACG...  
>NOT_PASS|CL000001|fdab42782e489ac806d35fed5867058b|V_16961|S_56|P_2|r_9.22178189031|SEED  
TACG...  
>p__Pseudomonadota|CL000005|6e105f3f79341d4c6ac024b928786898|V_1211301|S_8765|P_81|r_12.7105402863|SEED  
TACG...  

`makeblastdb -in SEEDS_97.0_WORKING_NAMES.fa -dbtype 'nucl' -out SEEDS_97.0_WORKING_NAMES`


### BINNING SINGLETONS

`zgrep '>' GB_VOL1_20250526_CLEANED_ranked_single.fa.gz | wc -l`
321663949

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`

`python2.7 split_fasta_by_group_size.py GB_VOL1_20250526_CLEANED_ranked_single.fa.gz 1290000`

BLAST SINGLETONS  

```
for file in *.fas
 do
  echo "blastn -query ${file} -db /mnt/DATA1/GLOBAL_BACTERIA/FINAL/BINNING/SEEDS_97.0_WORKING_NAMES -out ${file%%.fas}_SEEDS97.txt -outfmt 6 -evalue 1E-5 -num_threads 1 -max_target_seqs 10"
 done > blast_command.sh

cat blast_command.sh | parallel 
```

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db SEEDS_97.0_WORKING_NAMES -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```


### TAXONOMY FOR CLUSTERS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/get_taxonomy_for_clusters.py`

`python2.7 get_taxonomy_for_clusters.py GB_VOL1_20250526_CLEANED_uniq_Greengenes2_20250505_PROCESSED.txt.gz taxonomy.tsv SEEDS_97.0_WORKING_NAMES.fa SEEDS_97.0_WORKING_NAMES_TAXONOMY.txt`  


### SORT AND RENAME CLUSTERS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/sort_and_rename_perm_clusters.py`

process clusters and set version

`python2.7 sort_and_rename_perm_clusters.py SEEDS_97.0_WORKING_NAMES.fa 1`

### CREATE VARIANTS TABLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/create_variants_table_sample_pairs.py`

`python2.7 create_variants_table_sample_pairs.py GB_VOL1_20250526_CLEANED.fa.gz CLUSTERED_VARS_RANDOM_10000.fa SEEDS_97.0_WORKING_NAMES.fa.info`

OUTPUTS:  
VARIANTS_TABLE_CLUSTER_PAIRS.txt  
VARIANTS_TABLE_SAMPLE_PAIRS.txt  
VARIANTS_TABLE.txt.gz  

### CONVERT VARIANTS TABLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/convert_variants_table.py`

`python2.7 convert_variants_table.py VARIANTS_TABLE_CL.txt.gz`

### CREATE TAXONOMY TABLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/taxonomy_for_clusters.py`

`python2.7 taxonomy_for_clusters.py VARIANTS_TABLE_CLUSTER_PAIRS.txt SEEDS_97.0_WORKING_NAMES_TAXONOMY.txt`

OUTPUTS:  
TAXONOMY_CLUSTERS.txt  
TAXONOMY_SP_GEN.txt  

CLUSTERS - SPECIES - GENUS ABUNDANCE TABLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/create_taxonomy_tables.py`

`python2.7 create_taxonomy_tables.py VARIANTS_TABLE_TEST.txt TAXONOMY_CLUSTERS.txt`

OUTPUTS:  
ABUND_TABLE_CLUSTERS.txt  
ABUND_TABLE_SPECIES.txt  
ABUND_TABLE_GENERA.txt  

### SAMPLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/GENERATE_SAMPLES_FASTA_GB1.py`

`mkdir SAMPLES`

`python2.7 GENERATE_SAMPLES_FASTA_GB1.py VARIANTS_TABLE.txt.gz DATABASE_TABLES_NO_SINGLETONS/VARIANTS_TABLE_SAMPLE_PAIRS.txt SAMPLES/`

```
for file in *.fasta
 do zip -j ${file%%.fas}.zip $file
done

for file in *.fasta 
do
 rm -rf $file
done
```

### CREATE SAMPLE AND PAPER TABLES

count sample sequences from splited samples
```
for f in *.fasta.zip; do
   sum=$(unzip -p "$f" | grep "^>" | sed -E 's/.*size=([0-9]+).*/\1/' | awk '{s+=$1} END{print s+0}')
   echo -e "$f\t$sum"
 done > sample_counts.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/gb1_samples_metadata.txt.gz`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/VARIANTS_TABLE_SAMPLE_PAIRS.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/GB_PROCESS_SAMPLES.py`

`python2.7 GB_PROCESS_SAMPLES.py gb1_samples_metadata.txt.gz VARIANTS_TABLE_SAMPLE_PAIRS.txt`

OUTPUTS:  
SAMPLES_ADVANCED.txt
SAMPLES_BASIC.txt
SAMPLES_PAPERS.txt

### ALL FILES NEEDED FOR DATABASE

- on database server
`mkdir /mnt/data/mysql-data/GB1_TABLES_RAW`
`sudo chmod -R 777 /mnt/data/mysql-data/GB1_TABLES_RAW`

- copy files to database server
`scp -i /mnt/DATA1/KEYS/xxx.pem TAXONOMY_CLUSTERS.txt ubuntu@xxx.xxx.xxx.xxx:/mnt/data/mysql-data/GB1_TABLES_RAW`


#####################
### SERVER SET UP ###
#####################

### DOCKER

Install
`sudo apt update`
`sudo apt upgrade -y`
`sudo apt install -y ca-certificates curl gnupg lsb-release`
`sudo install -m 0755 -d /etc/apt/keyrings`
`curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg`
`sudo chmod a+r /etc/apt/keyrings/docker.gpg`
`echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] \
  https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" \
  | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null`
`sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin`

Check installation
`docker --version`
`sudo docker run hello-world`

If you want to see all containers (running + stopped):
`sudo docker ps -a`

### MariaDB - MySQL

# create a persistent volume
`docker volume create mariadb_data`

# run MariaDB with persistence and basic init
`docker run -d --name mariadb_ok \
  --restart unless-stopped \
  -p 3306:3306 \
  -e TZ=Europe/Prague \
  -e MARIADB_ROOT_PASSWORD='yourStrongRootPwd' \
  -e MARIADB_DATABASE='GB1' \
  -v mariadb_data:/var/lib/mysql \
  mariadb:12`

# watch startup
docker logs -f mariadb_ok

# RESTART MariaDB

`docker ps -a`
CONTAINER ID   IMAGE
2d19a8877979   mariadb:11

`docker container prune`  

`docker stop <container-id>`
`docker rm <container-id>`

# kill all dockers
`docker kill $(docker ps -aq)`
`docker rm mariadb_ok`

- be sure you mkdir `/mnt/data/mysql-data/GB1`

# RUN DATABASE
```
docker rm -f mariadb_ok 2>/dev/null

docker run \
  --name mariadb_ok \
  --mount type=bind,source=/mnt/data/mysql-data,target=/var/lib/mysql \
  -e MARIADB_ROOT_PASSWORD=root \
  -e MYSQL_USER=test \
  mariadb:12 \
  --innodb-buffer-pool-size=64G \
  --innodb-buffer-pool-instances=8 \
  --tmpdir=/var/lib/mysql/GB1
```

# Connect from host:
`docker exec -it mariadb_ok mariadb -u root -p`
[root]

```
use GB1;
SELECT * FROM  maillist LIMIT 20;
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
(1, 'GlobalBacteria', 'v1.0', '1', '2 (4.2.2020)', 51603682, 1116074647, 'First release', 'Větrovský T., Kyselková M., Baldrian P.: GlobalBacteria, a global database of bacterial occurrences from high-throughput-sequencing metabarcoding studies.', '5.8.2025');
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

```
CREATE TABLE IF NOT EXISTS `samples_basic` (
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
# TRUNCATE TABLE `samples_basic`;

```
LOAD DATA INFILE '/var/lib/mysql/GB1_TABLES_RAW/SAMPLES_BASIC.txt'
INTO TABLE `samples_basic`
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

ALTER TABLE samples_basic ADD INDEX(id);

ALTER TABLE samples_basic CHANGE ITS_total seqs_total INT;

```
CREATE TABLE IF NOT EXISTS `samples_advanced` (
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
LOAD DATA INFILE '/var/lib/mysql/GB1_TABLES_RAW/SAMPLES_ADVANCED.txt'
INTO TABLE `samples_advanced`
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

ALTER TABLE samples_advanced ADD INDEX(id);

```
CREATE TABLE IF NOT EXISTS `samples_papers` (
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

`LOAD DATA LOCAL INFILE '/var/lib/mysql/GB1_TABLES_RAW/SAMPLES_PAPERS.txt' INTO TABLE samples_papers FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

`ALTER TABLE samples_papers ADD INDEX(id);`

```
CREATE TABLE IF NOT EXISTS `clusters_tax` (
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

`LOAD DATA LOCAL INFILE '/var/lib/mysql/GB1_TABLES_RAW/TAXONOMY_CLUSTERS.txt' INTO TABLE clusters_tax FIELDS TERMINATED BY '\t' ESCAPED BY '\b';`

```
ALTER TABLE clusters_tax ADD INDEX(id);
ALTER TABLE clusters_tax ADD INDEX(cluster);
ALTER TABLE clusters_tax ADD INDEX(Species);
ALTER TABLE clusters_tax ADD INDEX(Genus);
CREATE INDEX idx_clusters_tax_species_id ON clusters_tax (Species, id);
CREATE INDEX idx_clusters_tax_genus_id ON clusters_tax (Genus, id);
ANALYZE TABLE clusters_tax, samplevar;
```

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

### SET APP USER

GRANT ALL privileges ON GB1.* TO 'test'@'%';
FLUSH PRIVILEGES;

SELECT User,Host FROM mysql.user;
REVOKE ALL PRIVILEGES, GRANT OPTION FROM 'test'@'%';
DROP USER 'test'@'%';
CREATE USER 'test' IDENTIFIED BY 'ubuntu';
GRANT ALL privileges ON GB1.* TO 'test'@'%';
FLUSH PRIVILEGES;

CREATE USER 'test'@'localhost' IDENTIFIED BY 'ubuntu';
GRANT ALL privileges ON GB1.* TO 'test'@'localhost';
FLUSH PRIVILEGES;

SHOW GRANTS FOR 'test'@'localhost';

### Docker app

/home/ubuntu/Docker_GlobalFungi

check Dockerfile

```
mv /home/ubuntu/FIX/* /home/ubuntu/Docker_GlobalFungi/app/

cd /home/ubuntu/Docker_GlobalFungi/

docker build --network host -t fungi_test . 
```

test app
`docker run --network host --mount type=bind,source=/mnt/data/databases_docker,target=/home/fungal/databases fungi_test R -e 'shiny::runApp(appDir = "/srv/shiny-server")'`

check it out
 `sensible-browser` 

Use a terminal browser

```
sudo apt-get update
sudo apt-get install -y w3m   # or: lynx, links2
w3m http://127.0.0.1:3838
```


### Docker VOLUME

mkdir /mnt/data/databases_docker

#docker volume create -d local-persist -o mountpoint=/mnt/data/databases_docker --name=databases
```
docker volume create \
  --driver local \
  --opt type=none \
  --opt o=bind \
  --opt device=/mnt/data/databases_docker \
  databases
```

`docker volume inspect databases`

docker volume inspect databases
[
    {
        "CreatedAt": "2025-08-21T23:10:32Z",
        "Driver": "local",
        "Labels": null,
        "Mountpoint": "/var/lib/docker/volumes/databases/_data",
        "Name": "databases",
        "Options": {
            "device": "/mnt/data/databases_docker",
            "o": "bind",
            "type": "none"
        },
        "Scope": "local"
    }
]

`docker container create --name temp -v databases:/home/fungal/databases busybox`

`docker run --volume databases:/home/fungal/databases -it fungi_test bash`
Leave and stop the container [Ctrl-D]

### SET THE NGINX

`nginx -t`

`/etc/nginx/nginx.conf`

### RUN THE DATABASE

stop testing proadcast
`sudo systemctl stop shiny-server`

sudo apt-get update
sudo apt-get install -y temurin-8-jdk
  
`java -jar shinyproxy-2.1.0.jar -Xloggc:shinyproxy_loggc.txt -Xmx12000m &> out_shinyproxy-2.1.0_20250000_12_31_01.txt`  
  
NEW!!!  
`java -Xms1g -Xmx12g -Xlog:gc*:file=shinyproxy_gc.log:time,uptime,level,tags -jar shinyproxy-3.2.0.jar`

### FIXING MySQL OVER

# 1) Restart DB (graceful)
docker restart mariadb_ok

# 2) Also restart the app container to reset its connection pool
docker restart crazy_jang

# Basic ping / processlist
docker exec -it mariadb_ok bash -lc 'mariadb-admin -uroot -p ping || true'  
docker exec -it mariadb_ok bash -lc 'mariadb -uroot -p -e "SHOW FULL PROCESSLIST;"'

Reduce “aborted connection” noise if clients are chatty/slow:
`SET GLOBAL wait_timeout=600, net_read_timeout=120, net_write_timeout=120;`
  
Warm it up  
```









