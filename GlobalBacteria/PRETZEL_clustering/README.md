### RANK AND SORT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/score_variants.py`

`python2.7 score_variants.py GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz`

Processed sequences 2001205289 to variants 514364177 (GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz)

### CLUSTERING

`/home/pretzel/ali`

```
#  /mnt/DATA1/Align/align512
IjaClust ver 3.2
Usage: align <[dir/]input.fasta[.gz|.gzip]> [similarity_threshold] [-o output_dir] [serial] [-m max_threads]
Example: align data/sequences.fasta.gz 97 -o output -m 8
Threshold can be defined with one decimal place. Default value is 97.0
Output file is written to input or provided directory with name based on input file name and threshold
Max threads can be set to 'serial' for single-threaded execution or a positive integer for multi-threading (default is 200 or number of hardware threads)
```

`/mnt/DATA1/Align/align512 GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.gz 97.0 -m 512`

IjaClust ver 3.2  
Path: GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.gz zipped Similarity threshold: 97 Max threads: 512  


### RESULTS

number of variants  514 364 177;  
created 87684198 clusters;  
duration 73.48h  
(total sequences  2 001 205 289)  

```
(
echo -e "clusterName\tnumberOfVariants\tnumberOfSequences"
zgrep '^>' GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz \
| awk -F'[|;=]' '
{
    cluster = substr($1, 2)
    size = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "size") {
            size = $(i+1)
            break
        }
    }
    count[cluster]++
    sum[cluster] += size
}
END {
    for (c in count) {
        print c "\t" count[c] "\t" sum[c]
    }
}' | sort
) > cluster_stats.tsv
```

number of singletons  
`awk 'NR>1 && $2==1 && $3==1 {n++} END {print n+0}' cluster_stats.tsv`

singletons (clusters with one variant of size 1) 69 754 077  
17 930 121 clusters with at least doubletons  

`awk -F '\t' 'NR>1 {print $3}' cluster_stats.tsv | sort | uniq -c > cluster_size_info.txt`

### GET OTU TAB - REDUCED FORMATE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/GET_OTUTAB_REDUCED_FORMATE.py`

```
python2.7 GET_OTUTAB_REDUCED_FORMATE.py GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz
```

### IDENTIFICATION OF CLUSTERS

`# wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/select_cluster_representatives.py`

```
# python select_cluster_representatives.py \
#   -i GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz \
#   -o cluster_representatives.fa \
#   --seed 123
```

get SEED sequence from cluster

```
zcat GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz \
  | tr -d '\r' \
  | grep -A 1 --no-group-separator '|100\.0$' > GB_BOTH_VOL_20260413_97_clustered_SEEDs.fa
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`


```
python2.7 split_fasta_by_group_size.py GB_BOTH_VOL_20260413_97_clustered_SEEDs.fa 350800
```


`makeblastdb -in Greengenes2_2024_09_backbone_taxonomy.fas -dbtype 'nucl' -out greenegenes2_2024_09`

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db greenegenes2_2024_09 -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```

```
awk 'BEGIN{FS=OFS="\t"} {sub(/;.*/, "", $1); print}' GB_BOTH_VOL_20260413_97_clustered_SEEDs_MISSED_COMPLETE_best.txt > SEEDs_BLAST_COMPLETE_greenegenes2_2024_09_best.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/PROCESS_CLUSTER_REP_BLAST_NEW.py`

```
python2.7 PROCESS_CLUSTER_REP_BLAST_NEW.py GB_BOTH_VOL_20260413_97_clustered_SEEDs.fa.gz SEEDs_BLAST_COMPLETE_greenegenes2_2024_09_best.txt SEEDs_BLAST_COMPLETE_greenegenes2_2024_09_PROCESSED.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/get_cluster_info_final.py`

```
python2.7 get_cluster_info_final.py GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz CLUSTERS_IDENT_greenegenes2_2024_09_PROCESSED.txt GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz samples_and_studies.txt CLUSTERS_INFO_FINAL.txt
```
### SILVA

`makeblastdb -in silva_138_1_fixed_taxa_MithoChloro_FINAL_20230818.fasta -dbtype 'nucl' -out SILVA_138_1_MitoChloro`

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas
do
 echo "blastn -query ${file} -db SILVA_138_1_MitoChloro -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas}_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```

### VARIANTS TABLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/create_variants_table_and_sample_pairs.py`

```
python2.7 create_variants_table_and_sample_pairs.py GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz
```

### TAXONOMY TABLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/taxonomy_for_clusters_new.py`
  
https://ftp.microbio.me/greengenes_release/current/
```
python2.7 taxonomy_for_clusters_new.py CLUSTERS_INFO_TOTAL_FINAL.txt.gz 2024.09.taxonomy.id.tsv.gz 2
```

### CLUSTERS - SPECIES - GENUS ABUNDANCE TABLES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/create_taxonomy_tables.py`

`python2.7 create_taxonomy_tables.py VARIANTS_TABLE.txt.gz TAXONOMY_CLUSTERS.txt`

OUTPUTS:  
ABUND_TABLE_CLUSTERS.txt  
ABUND_TABLE_SPECIES.txt  
ABUND_TABLE_GENERA.txt  

`awk 'BEGIN{OFS="\t"} {$1=$1".2"; print}' ABUND_TABLE_CLUSTERS.txt > ABUND_TABLE_CLUSTERS_v2.txt`

### UPDATE TAXONOMY CLUSTERS

```
gzip -dc TAXONOMY_CLUSTERS.txt.gz |
awk 'BEGIN { OFS="\t" } {
  $1 = $1 ".2"
  print NR, $0
}' |
gzip > TAXONOMY_CLUSTERS_GOOOD.txt.gz
```

### SIMPLIFY CLUSTER NAMES

`zcat VARIANTS_TABLE.txt.gz | awk 'BEGIN{OFS="\t"} {sub(/^CL0*/, "", $4); print}' | gzip > VARIANTS_TABLE_CLNUM.txt.gz`

### VARIANTS FASTA

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/GB_GET_VARINATS_FOR_ALL_TAXA_FINAL.py`

`mkdir VARIANTS`

`python2.7 GB_GET_VARINATS_FOR_ALL_TAXA_FINAL.py VARIANTS_TABLE_CLNUM.txt.gz TAXONOMY_CLUSTERS_GOOOD.txt.gz VARIANTS/`

### CONVERT VARIANTS TABLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PermanentClusters/convert_variants_table.py`

`python2.7 convert_variants_table.py VARIANTS_TABLE_CLNUM.txt.gz`

VARIANTS_variants.txt.gz - columns: id(int); cl_id(int); hash(varchar); sequence(text);  
VARIANTS_samplevar.txt.gz - columns: id(bigint); variant(int); sample(int); abundance(int); cl_id(int);   


### BINNING NEW SEQUENCES

sample names at the beggining  
```
gzip -dc Bacteria_joined_qm30_trimmed_noambi.fa.gz | awk '/^>/ {
  split($0, a, "|")
  print ">" a[2] "|" substr(a[1], 2)
  next
}
{ print }' | gzip > Bacteria_joined_qm30_trimmed_noambi_samplefirst.fa.gz
```

generate scored variants   
```
python2.7 score_variants.py Bacteria_joined_qm30_trimmed_noambi_samplefirst.fa.gz
```

combine the original seed representatives with scored variants  
```
cat /mnt/DATA/projects/avetrot/GLOBAL_BACTERIA_BOTH_VOLUMES/FINAL/FINAL_TABLES/GTDB/GB_VOL2_AS_DATABASE/GB_BOTH_VOL_20260413_97_clustered_SEEDs.short.fa.gz  Bacteria_joined_qm30_trimmed_noambi_samplefirst.fa.gz_scored_variants.fa.gz > Bacteria_scored_variants_for97sim_clusters.fa.gz
```
cluster all together  
```
/mnt/DATA1/Align/align512 Bacteria_scored_variants_for97sim_clusters.fa.gz 97.0 -m 512
```
remove the original seed representatives  
```
gzip -dc Bacteria_scored_variants_for97sim_clusters.fa.97.clustered.gz | awk '
  /^>/ {
    keep = ($0 !~ /\|CL/)
  }
  keep
' | gzip > Bacteria_scored_variants_for97sim_clusters.fa.97.clustered_onlyDeadWood.gz
```
make reduced otu table  
```
python make_otu_table_multi.py \
  --cl-vars GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz FOR_SAILEE/Bacteria_scored_variants_for97sim_clusters.fa.97.clustered_onlyDeadWood.gz \
  --fasta GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz FOR_SAILEE/Bacteria_joined_qm30_trimmed_noambi_samplefirst.fa.gz \
  -o FOR_SAILEE/compressed_otu_tab_otus_to_samples.txt
```

### GTDB

binning to extracted ssu  
```
/mnt/DATA1/Align/align512 ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.gz 97.0 -m 320
```
Total clusters: 87694333  
Timings: read=40.85min cluster=51.74h total=52.42h  
cluster file: ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered.gz  

remove the original seed representatives 
```
gzip -dc ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered.gz | awk '
>   /^>/ {
>     keep = ($0 !~ /\|CL/)
>   }
>   keep
> ' | gzip > ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered_onlyGTDB.gz
```
get only GTDB ssu binned to existing GB clusters
```
gzip -cd ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered_onlyGTDB.gz |
awk -v max=87684198 '
BEGIN { RS=">"; ORS="" }
NR > 1 {
    header = substr($0, 1, index($0, "\n") - 1)
    cluster = header
    sub(/^CL/, "", cluster)
    sub(/\|.*/, "", cluster)

    if ((cluster + 0) <= max)
        print ">" $0
}
' |
gzip > ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered_onlyGTDB_onlyGBclusters.gz
```

GTDB 97sim for clusters  
```
zgrep '>' ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered_onlyGTDB_onlyGBclusters.gz | awk -F'>' '{print $2}' | awk -F'~' '{print $1}' | awk -F'|' '{print $1"\t"$2}' | sort | uniq > FINAL/RAW_GTDB_CLUSTERS97sim.txt
```
GTDB to cluster variants representatives  
```
zcat ssu_all_r232_GB_EXTRACTED_for97sim_clusters.fa.97.clustered_onlyGTDB_onlyGBclusters.gz \
  | tr -d '\r' \
  | grep --no-group-separator '|100\.0$' | awk -F'>' '{print $2}' | awk -F'~' '{print $1}' | awk -F'|' '{print $1"\t"$2}' | sort | uniq > FINAL/RAW_GTDB_CLUSTERS_vars.txt
```

```
wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/dereplicate_to_md5_gz.py
python2.7 dereplicate_to_md5_gz.py ssu_all_r232_GB_EXTRACTED.fas.gz ssu_all_r232_GB_EXTRACTED_md5_derep.fas ssu_all_r232_GB_EXTRACTED_md5_derep.tab
grep '>' ssu_all_r232_GB_EXTRACTED_md5_derep.fas | awk -F'>' '{print $2}' > FINAL/GTDB_md5.txt
zgrep --no-group-separator -F -f GTDB/FINAL/GTDB_md5.txt VARIANTS_variants.txt.gz > GTDB/FINAL/VARIANTS_variants_GTDBonly.txt
awk -F'~' '{print $1}' ssu_all_r232_GB_EXTRACTED_md5_derep.tab | sort | uniq > FINAL/ssu_all_r232_GB_EXTRACTED_md5_uniq_acc.txt
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/get_gtdb_for_exact_vars.py`
`python2.7 get_gtdb_for_exact_vars.py VARIANTS_variants_GTDBonly.txt ssu_all_r232_GB_EXTRACTED_md5_uniq_acc.txt GB2_GTDB_variants.txt`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/PRETZEL_clustering/get_gtdb_for_clusters_97sim.py`
`python2.7 get_gtdb_for_clusters_97sim.py RAW_GTDB_CLUSTERS97sim.txt GB2_GTDB_CLUSTERS97sim.txt`







