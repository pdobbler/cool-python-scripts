# Amplicon data

## GET THE DATA

`wget -O AMPLICON_DATA.zip "http://www.biomed.cas.cz/mbu/lbwrf/example/AMPLICON_DATA.zip"`

```
unzip AMPLICON_DATA.zip
cd AMPLICON_DATA/
```

## OTU BASED - CLUSTERING APPROACH

`cd 16S_example_data/`

**PREPARATION OF SAMPLE TABLE FROM FWD and REV tags**
```
head samples.txt
head barcoded_primers_fwd.txt
head barcoded_primers_rev.txt
```

**be sure about formate**
```
dos2unix samples.txt
dos2unix barcoded_primers_fwd.txt
dos2unix barcoded_primers_rev.txt
```

```
awk -F'\t' '{print $2}' samples.txt > fwd_tags.txt
awk -F'\t' '{print $3}' samples.txt > rev_tags.txt
```


### CREATE AWK SCRIPT

`nano vlookup.awk`

```
# vlookup script
FNR==NR{
  a[$1]=$col
  next
}
{OFS="\t"; if ($1 in a) {print a[$1], $1} else {print "NA", $1}  }
```

**Ctrl+o & Enter (save)**

**Ctrl+x         (exit)**



**run the script**
```
awk -v col='3' -f vlookup.awk barcoded_primers_fwd.txt fwd_tags.txt > cols_fwd.txt
awk -v col='3' -f vlookup.awk barcoded_primers_rev.txt rev_tags.txt > cols_rev.txt
awk -F'\t' '{print $1}' samples.txt > cols_sample.txt
```

**function "paste" joins the columns of the table**

```
paste -d'\t' cols_sample.txt cols_fwd.txt cols_rev.txt > barcode_table.txt
```

**or**

```
paste <(awk -F'\t' '{print $1}' samples.txt) <(awk -v col='3' -f vlookup.awk barcoded_primers_fwd.txt fwd_tags.txt) <(awk -v col='3' -f vlookup.awk barcoded_primers_rev.txt rev_tags.txt)
```


### DEMULTIPLEXING

**move data to a new folder**

```
mkdir demulti
mv barcode_table.txt demulti/
mv *.fastq demulti/
cd demulti/
```

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Amplicons/Course_MB140P94/demultiplexing_fastq.py`

`python2.7 demultiplexing_fastq.py barcode_table.txt BAC_R1.fastq BAC_R2.fastq`

**count sample sequences**

`grep '^@M' BAC_R1.fastq.cut | awk -F'|' '{print $2}' | sort | uniq -c > bacteria_counts.txt`


### JOINNING PAIRED-ENDS - fastq-join

**move data to a new folder**

```
mkdir joining
mv *.cut joining/
cd joining/
```

`fastq-join -v " " -p 15 -m 40 BAC_R1.fastq.cut BAC_R2.fastq.cut -o BAC_joined`


### QUALITY FILTERING

**average quality qreater or equal than treshold**

`seqkit seq -Q 30  BAC_joinedjoin > BAC_joined_qm30.fq`

**convert FASTQ to FASTA**

`seqkit fq2fa BAC_joined_qm30.fq > BAC_joined_qm30.fa`


### TRIMM Bacterial 16S - rest of primers & lenght restriction 

```
mkdir trimm
mv BAC_joined_qm30.fa trimm/
cd trimm/
```

**nucleotides only - Max**

`awk 'NR % 2 == 0' BAC_joined_qm30.fa | wc -L`

**nucleotides only - Min**

`awk 'NR % 2 == 0' BAC_joined_qm30.fa | awk '{print length}' | sort -n | head -n1`

**average length**

```
chars=$(awk 'NR % 2 == 0' BAC_joined_qm30.fa | wc -c)
words=$(awk 'NR % 2 == 0' BAC_joined_qm30.fa | wc -w)
avg_word_size=$(( ${chars} / ${words} ))
echo "average length $avg_word_size"
```

**median**

`sort -n <(awk 'NR % 2 == 0' BAC_joined_qm30.fa | awk '{ print length }') |awk '{a[NR]=$0}END{print(NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}'`


### trimm primer residuals

**16S - primers: 515F & 806R**
**284 - (15 + 16) = 253 +- 8 bp > length range**


**remove first 15 bases - forward primer residual**

`seqkit subseq -r 16:-1 BAC_joined_qm30.fa > BAC_joined_qm30_f.fa`

**remove last 16 bases - reverse primer residual**

`seqkit subseq -r 1:-17 BAC_joined_qm30_f.fa > BAC_joined_qm30_fr.fa`

**sequences shorter than or equal to the maximum length**

`seqkit seq -g -m 245 BAC_joined_qm30_fr.fa > BAC_joined_qm30_fr_min245.fa`

**sequences longer than or equal to the minimum length**

`seqkit seq -g -M 261 BAC_joined_qm30_fr_min245.fa > BAC_joined_qm30_fr_min245_max261.fa`

**linearize**

`seqkit seq -w 0 BAC_joined_qm30_fr_min245_max261.fa > BAC_joined_qm30_trimmed_linear.fa`


**length Max & Min**

```
awk 'NR % 2 == 0' BAC_joined_qm30_trimmed_linear.fa | wc -L
awk 'NR % 2 == 0' BAC_joined_qm30_trimmed_linear.fa | awk '{print length}' | sort -n | head -n1
```

