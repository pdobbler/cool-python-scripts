### Annotate new variants

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/dereplicate_to_md5_gz.py`

`python2.7 dereplicate_to_md5_gz.py itsx_REL6_ITSX_all_parts_ITS1_final.fa.gz itsx_REL6_ITSX_all_parts_ITS1_md5_variants.fa itsx_REL6_ITSX_all_parts_ITS1_derep.map`


580970909 sequences loaded correctly - 0 sequnces are empty - Omitted!
580970909 sequence variants sorted ...
Dereplication is done - 76579647 groups from 580970909 seqs

### ANNOTATE VARIANTS

[root@hpe ITSx]# zgrep '>' GF6_ALL_SAMPLES.fa.its1.gz_scored_variants.fa.gz | wc -l
274525757

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`


```
python2.7 split_fasta_by_group_size.py GF6_ALL_SAMPLES.fa.its1.gz_scored_variants.fa.gz 1100000
```

[root@hpe ITSx]# zgrep '>' GF6_ALL_SAMPLES.fa.its2.gz_scored_variants.fa.gz | wc -l
490801868

```
python2.7 split_fasta_by_group_size.py GF6_ALL_SAMPLES.fa.its2.gz_scored_variants.fa.gz 2000000
```

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db /mnt/DATA/DATABASES/UNITE10/UNITE_10_SIMPLE -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```
  
### get best hits
  
`cat ANNOT_ITS1_UNITE10/*_best.tab > GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_best.tab`
`cat ANNOT_ITS2_UNITE10/*_best.tab > GF6_ALL_SAMPLES_its2_scored_variants_UNITE10_best.tab`

  
### get processed blast file

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/PROCESS_BLAST_RESULTS_SIMPLE.py`

`python2.7 PROCESS_BLAST_RESULTS_SIMPLE.py GF6_ALL_SAMPLES.fa.its1.gz_scored_variants.fa.gz GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_best.tab GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_PROCESSED.txt`

```
awk 'BEGIN{OFS="\t"} NR==1{print; next} {sub(/;.*/, "", $1); print}' \
GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_PROCESSED.txt \
> GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_PROCESSED_clean.txt
```
### CREATE VARIANTS TABLE
  
`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/RELESE6/create_variant_table.py`
  
```
python3 create_variant_table.py \
    --fasta GF6_ALL_SAMPLES.renamed.fa.its1.gz \
    --table GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_PROCESSED_clean.txt.gz \
    --similarity 98.5 \
    --coverage 90.0 \
    --marker ITS1 \
    --output GF6_ALL_SAMPLES_ITS1_variants.tsv
```

```
Reading FASTA: GF6_ALL_SAMPLES.fa.its1.gz
FASTA records: 1845276525
Unique sequence variants: 274525757
Reading identification table: GF6_ALL_SAMPLES_its1_scored_variants_UNITE10_PROCESSED_clean.txt.gz
Identification table rows: 274525757
Unique QUERY IDs in identification table: 274525757
QUERY IDs without HIT: 50605033
Writing output: GF6_ALL_SAMPLES_ITS1_variants.tsv

=== FINAL REPORT ===
FASTA records:                  1845276525
Unique sequence variants:       274525757
Identification table rows:      274525757
Unique identification QUERYs:   274525757
Passed thresholds:              27611454
Failed thresholds:              196309270
QUERYs without HIT (-):         50605033
FASTA seqIDs without QUERY:     0
Duplicated QUERY IDs:           0

No missing or duplicated QUERY IDs detected.

Finished.
```

```
Reading FASTA: GF6_ALL_SAMPLES.fa.its2.gz
FASTA records: 4579522152
Unique sequence variants: 490801868
Reading identification table: GF6_ALL_SAMPLES_its2_scored_variants_UNITE10_PROCESSED_clean.txt.gz
Identification table rows: 490801868
Unique QUERY IDs in identification table: 490801868
QUERY IDs without HIT: 34835151
Writing output: GF6_ALL_SAMPLES_ITS2_variants.tsv

=== FINAL REPORT ===
FASTA records:                  4579522152
Unique sequence variants:       490801868
Identification table rows:      490801868
Unique identification QUERYs:   490801868
Passed thresholds:              65904385
Failed thresholds:              390062332
QUERYs without HIT (-):         34835151
FASTA seqIDs without QUERY:     0
Duplicated QUERY IDs:           0

No missing or duplicated QUERY IDs detected.

Finished.
```
  
### DATABASE STRUCTURE
  
```
CREATE TABLE IF NOT EXISTS `variants` (
  `id` int(10) unsigned NOT NULL,
  `cl_id` int(10) unsigned NOT NULL,
  `marker` varchar(4) NOT NULL,  
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
