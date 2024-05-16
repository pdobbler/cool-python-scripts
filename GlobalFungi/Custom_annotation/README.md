
### split variants into groups

`python2.7 split_fasta_by_group_size.py GF5_RAW_TABLE_PROCESSED_VARIANTS.fa 2400000`


### make custom annotation database

`makeblastdb -in ErMF_isolates_seq_NCBI.fas -parse_seqids -dbtype nucl`


### blast

```
for file in *.fas
do  
 echo "blastn -query ${file} -db ErMF_isolates_seq_NCBI -out ${file%%.fas}.ErMF_NCBI.txt -evalue 1E-5 -outfmt 6 -num_threads 1 -max_target_seqs 10"
done > blast_command.sh
```

`cat blast_command.sh | parallel`


### sort

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

```
for file in *.txt
do  
 echo "sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr ${file} | sort -u -k1,1 --merge > ${file%%.txt}_best.tab"
done > sorting.sh
```

`cat sorting.sh | parallel`


### combine best results

`cat *_best.tab > GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt`


### extract best hit variants


`awk -F'\t' '{print $1}' GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt > GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best_SEQNAMES.txt`

`zgrep --no-group-separator -A 1 -F -f GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best_SEQNAMES.txt GF5_RAW_TABLE_PROCESSED_VARIANTS.fa.gz > ErMF_NCBI/GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI.fa`

### get processed blast file

`python2.7 PROCESS_BLAST_RESULTS_SIMPLE.py GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI.fa GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_PROCESSED.txt`


### generate otu-table


`python2.7 GET_TABLE_FROM_PROCESSEDBL_AND_VARSTAB.py GF5_RAW_TABLE_PROCESSED.txt.gz ErMF_NCBI/GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_PROCESSED.txt ErMF_NCBI/ErMF_NCBI_ABUND_TABLE.txt 98.5 98 GF5_RAW_TABLE_SAMPLES.txt.gz`

