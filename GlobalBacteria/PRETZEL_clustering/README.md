### RANK AND SORT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PRETZEL_clustering/score_variants.py`

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


### IDENTIFICATION OF CLUSTERS

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PRETZEL_clustering/select_cluster_representatives.py`

```
python select_cluster_representatives.py \
  -i GB_BOTH_VOL_20260413_RENAMED_filtered.fa.gz_scored_variants.fa.97.clustered.gz \
  -o cluster_representatives.fa \
  --seed 123
```

get most abundant sequence from cluster


`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`


```
python2.7 split_fasta_by_group_size.py cluster_representatives.fa 350800
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

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/PROCESS_BLAST_RESULT.py`


```
python2.7 PROCESS_BLAST_RESULT.py cluster_representatives.fa cluster_representatives_greenegenes2_2024_09_best.txt cluster_representatives_greenegenes2_2024_09_PROCESSED.txt ITS2
```
