# Metagenome data

## GET THE DATA

`wget -O MG_DATA.tar.gz "http://www.biomed.cas.cz/mbu/lbwrf/example/mg_samples.tar.gz"`

```
tar -xvzf MG_DATA.tar.gz
rm MG_DATA.tar.gz
cd mg_samples/
```

## INTERLEAVE

### [KHMER](https://khmer.readthedocs.io/en/latest/) 
`conda activate khmerEnv`

```
for file in *_R1_001.fastq
do
   sample=${file%%_R1_001.fastq}
   echo "interleave-reads.py ${sample}_R1_001.fastq ${sample}_R2_001.fastq -o ${sample}.pe.fq"
done > interleave.sh
```
`cat interleave.sh | parallel`

### remove original files
rm *.fastq

### get stats

`seqkit stats *.pe.fq > stats_raw_seqs.txt`

`head stats_raw_seqs.txt`

## REMOVE ADAPTERS & SHORT READS

### [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

`conda activate cutadaptenv`

`cutadapt --version`

### [Index 1 (i7) Adapters](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/CDIndexes.htm)
### GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
### ATCTCGTATGCCGTCTTCTGCTTG

Check if they are there...

`grep 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' Sample01.pe.fq`

`grep 'ATCTCGTATGCCGTCTTCTGCTTG' Sample01.pe.fq`

```
for file in *.pe.fq
do
   sample=${file%%.pe.fq}
   echo "cutadapt -m 50 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g ATCTCGTATGCCGTCTTCTGCTTG -o ${sample}.pe.tr.fq ${file}"
done > trimming.sh
```

`cat trimming.sh | parallel &> out_trimming.txt`
           
### get stats

`seqkit stats *.pe.tr.fq > stats_trimmed_seqs.txt`

`head stats_trimmed_seqs.txt`


## QUALITY FILTERING

### [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
`conda activate fastx-env`

### Inspect read quality profiles before quality filtering
`fastx_quality_stats -Q33 -i Sample01.pe.tr.fq -o Sample01_raw.qstats.txt`
`fastq_quality_boxplot_graph.sh -i Sample01_raw.qstats.txt -o Sample01_raw.qstats.png`

```
for file in *.tr.fq
do
  sample=${file%%.tr.fq}		
  echo "fastq_quality_filter -i ${file} -Q33 -q 30 -p 85 -o ${sample}.tr.qc.fq"
done > qual_filter.sh
```

`cat qual_filter.sh | parallel &> out_qual.txt`

### Inspect read quality profiles after quality filtering

`fastx_quality_stats -Q33 -i Sample01.pe.tr.qc.fq -o Sample01_qc.qstats.txt`

`fastq_quality_boxplot_graph.sh -i Sample01_qc.qstats.txt -o Sample01_qc.qstats.png`

## remove short sequences

`mkdir filtered`

`mv *.qc.fq filtered`

`cd filtered/` 

### get stats

`seqkit stats *.qc.fq > stats_qc_seqs.txt`

`head stats_qc_seqs.txt`

```
for file in *.qc.fq
do
 sample=${file%%.qc.fq}
 echo "seqkit seq -g -m 250 ${file} > ${sample}.qc.cut"
done > cut_short.sh
```

`cat cut_short.sh | parallel &> out_cut.txt`

### get stats

`seqkit stats *.qc.cut > stats_cut_seqs.txt`

`head stats_cut_seqs.txt`

## extracting paired ends from the interleaved files

### activate khmer environment

`conda activate khmerEnv`

```
for file in *.qc.cut
do
   echo "extract-paired-reads.py ${file}"
done > extract_command.sh
```

`cat extract_command.sh | parallel`


`mkdir mapping`

`mv *.qc.cut mapping`

### mapping - these reads will be mapped

## rename files and merging se file

`mkdir renaming`

`mv *.qc.cut.* renaming`

`cd renaming`


```
for file in *.pe
do
   sample=${file%%.pe.tr.qc.cut.pe}
   mv ${file} ${sample}.pe.fq
done
```

```
for file in *.se
do
   sample=${file%%.pe.tr.qc.cut.se}
   mv ${file} ${sample}.se.fq
done
```
