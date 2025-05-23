
### split variants into groups

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/split_fasta_by_group_size.py`

`python2.7 split_fasta_by_group_size.py GF5_RAW_TABLE_PROCESSED_VARIANTS.fa 2400000`


### make custom annotation database

`makeblastdb -in ErMF_isolates_seq_NCBI.fas -parse_seqids -dbtype nucl`


### blast & sort

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.fas.gz
do
 echo "gunzip -c ${file} | blastn -query - -db OMF_sequences.fa -outfmt 6 -evalue 1E-5 -num_threads 2 -max_target_seqs 10 | sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge > ${file%%.fas.gz}_OMF_best.tab"
done > blast_and_sort_command.sh

mkdir -p /mnt/DATA1/tmp
export TMPDIR=/mnt/DATA1/tmp
cat blast_and_sort_command.sh | parallel --tmpdir /mnt/DATA1/tmp
```

### combine best results

`cat *_best.tab > GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt`


### extract best hit variants


`awk -F'\t' '{print $1}' GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt > GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best_SEQNAMES.txt`

`zgrep --no-group-separator -A 1 -F -f GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best_SEQNAMES.txt GF5_RAW_TABLE_PROCESSED_VARIANTS.fa.gz > ErMF_NCBI/GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI.fa`

### get processed blast file

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/PROCESS_BLAST_RESULTS_SIMPLE.py`

`python2.7 PROCESS_BLAST_RESULTS_SIMPLE.py GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI.fa GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_best.txt GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_PROCESSED.txt`


### generate otu-table

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/GET_TABLE_FROM_PROCESSEDBL_AND_VARSTAB.py`

`python2.7 GET_TABLE_FROM_PROCESSEDBL_AND_VARSTAB.py GF5_RAW_TABLE_PROCESSED.txt.gz ErMF_NCBI/GF5_RAW_TABLE_PROCESSED_VARIANTS_ErMF_NCBI_PROCESSED.txt ErMF_NCBI/ErMF_NCBI_ABUND_TABLE.txt 98.5 98.0 GF5_RAW_TABLE_SAMPLES.txt.gz`



## SPECIAL CASE WHEN ANOTATED DEREP SEQUENCES FROM SEED INSTEAD OF PROCESSED VARIANTS 
### generate otu-table from mapping file

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/Custom_annotation/GET_TABLE_FROM_PROCESSEDBL_AND_DEREPFILES.py`

`python2.7 GET_TABLE_FROM_PROCESSEDBL_AND_DEREPFILES.py NAKI_castleparks_vlk_corrected_derep.map NAKI_castleparks_clean_UNITE10_PROCESSED.txt NAKI_OTU_TAB_s98_5_c95.txt 98.5 95.0`

