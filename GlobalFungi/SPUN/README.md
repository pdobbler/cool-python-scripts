
### JOINING - PARALLEL WITH PERCENTAGE OF JOUNED READS IN joining.txt

```
output_file="joining.txt"
> "$output_file"

export output_file

parallel --bar --keep-order '
  sample={1}
  file1="${sample}_1.fastq"
  file2="${sample}_2.fastq"
  line_count=$(fastq-join -v " " -p 15 -m 40 "$file1" "$file2" -o "${sample}_joined" | awk "{ORS=\";\"; print} END {print \"\"}")

  total_reads=$(echo "$line_count" | grep -oP "Total reads: \K[0-9]+")
  total_joined=$(echo "$line_count" | grep -oP "Total joined: \K[0-9]+")

  if [[ -n "$total_reads" && -n "$total_joined" && "$total_reads" -gt 0 ]]; then
    percentage=$(awk "BEGIN {printf \"%.2f\", ($total_joined / $total_reads) * 100}")
  else
    percentage="N/A"
  fi

  echo "$sample ($file1, $file2);$line_count $percentage%" >> "$output_file"
' ::: $(ls *_1.fastq | sed 's/_1.fastq//')
```

### QUALITY FILTERING

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalBacteria/DATA_PROCESSING/filter_fastq_by_quality_mean_GZIP.py`

```
for file in *_joinedjoin; do
  sample=${file%%_joinedjoin};  
  echo "python2.7 filter_fastq_by_quality_mean_GZIP.py ${file} QM20/${sample}_qm20.fa 20 True"; 
done > filter_fastq.sh

cat filter_fastq.sh | parallel
```

### SELECT NON-EMPTY FILES

```
mkdir PROCESSED

for file in *.fa.gz; do
  if zgrep -q '^>' "$file"; then     mv "$file" PROCESSED/;   
  fi;
done
```

### ADD SAMPLE NAMES TO THE HEADERS


```
for file in *_qm20.fa.gz; do
    prefix=$(echo "$file" | cut -d'_' -f1,2)
    gzip -dc "$file" |     awk -v p="$prefix" '{ if ($0 ~ /^>/) { print ">" p "|" substr($0,2) } else { print $0 } }' |     gzip > "${file%.fa.gz}_modified.fa.gz"; 
done
```

### UPDATE SAMPLE NAMES AND PROJECT

original headers:  
>101_7069|VH01408:65:AAGHHTCM5:1:1101:58791:1303


`zcat SPUN2_joined_qm20_renamed.fa.gz |   sed 's/^>\(101_[0-9]*\)|/\>SP_\1|SPUN2|/' |   gzip > SPUN2_modified.fa.gz`

modified headers:  
>SP_101_7069|SPUN2|VH01408:65:AAGHHTCM5:1:1101:58791:1303


### ITSx

DEREPLICATE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Amplicons/Course_MB140P94/dereplicate_FASTA.py`

`python2.7 dereplicate_FASTA.py SPUN2_joined_qm20_renamed.fa.gz SPUN2_joined_qm20_renamed_derep.fasta SPUN2_joined_qm20_renamed_mapping.table`

`python2.7 dereplicate_FASTA.py SPUN2_modified.fa.gz SPUN2_joined_qm20_renamed_derep.fasta SPUN2_joined_qm20_renamed_mapping.table`

SPLIT

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/PermanentClusters/split_fasta_by_group_size.py`

`python2.7 split_fasta_by_group_size.py SPUN2_joined_qm20_renamed_derep.fasta.gz 550000`

EXTRACT ITS

`mv *.fas ITSx/`

```
for file in *.fas
do
 echo "/home/kdanielmorais/bioinformatics/tools/ITSx_1.0.11/ITSx -i ${file} --cpu 2 --only_full T -t F -o itsx_${file%%.fas}"
done > itsx_command.sh

cat itsx_command.sh | parallel &> out.txt
```

`cat ITSx/*.ITS2.fasta > SPUN2_joined_qm20_renamed_derep_ITS2.fasta`

RE-REPLICATE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Amplicons/Course_MB140P94/rereplicate_FASTA.py`

`python2.7 rereplicate_FASTA.py SPUN2_joined_qm20_renamed_derep_ITS2.fasta SPUN2_joined_qm20_renamed_mapping.table SPUN2_joined_qm20_renamed_ITS2.fasta`

### PREPARE CLUSTERING

zgrep '>' GF5_ALL_SAMPLES_AND_SPUN2_ITS2.fa.gz | wc -l
4 228 585 157


### MAKE VARIANTS WITH SIZE INFO and COMPARE WITH OLD ONES

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/clustering_usearch/make_uniques.py`

`python2.7 make_uniques.py SPUN2_joined_qm20_renamed_ITS2.fasta` 

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/GlobalFungi/SPUN/get_undetected_variants.py`








