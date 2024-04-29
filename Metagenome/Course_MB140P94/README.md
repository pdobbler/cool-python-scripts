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
<a href="https://www.openai.com" target="_blank">OpenAI</a>

### GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
### ATCTCGTATGCCGTCTTCTGCTTG

`grep 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' Sample01.pe.fq`
`grep 'ATCTCGTATGCCGTCTTCTGCTTG' Sample01.pe.fq`

