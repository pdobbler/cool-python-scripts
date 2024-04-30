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


## prepare for assembly

```
for file in *.pe.fq
do
   echo "split-paired-reads.py ${file}"
done > split_command.sh
```

`cat split_command.sh | parallel`

`mkdir forassembly`

`cat *.1 > forassembly/all.pe.fq.1`

`cat *.2 > forassembly/all.pe.fq.2`

`cat *.se.fq > forassembly/all.se.fq`


## assembly - [MEGAHIT](https://github.com/voutcn/megahit)
         
`cd forassembly`

`megahit -m 0.75 -t 8 -1 all.pe.fq.1 -2 all.pe.fq.2 -r all.se.fq -o all.Megahit.assembly`


## evaluating an assembly - [quast](https://github.com/ablab/quast)

`cd all.Megahit.assembly`

`quast  final.contigs.fa -o report`

### Check the report folder


## mapping reads to assembly

### making reference database

`cp final.contigs.fa ~/mg_samples/filtered/mapping/`
`cd  ~/mg_samples/filtered/mapping/`


`mkdir build`

`bowtie2-build final.contigs.fa build/final.contigs.build`

### aligning

```
for file in *.cut
do
 sample=${file%%.pe.tr.qc.cut}
 echo "bowtie2 -p 8 -x build/final.contigs.build -q ${file} -S ${sample}.sam"
done > align.sh
```

`cat align.sh | parallel`

### sam to bam

```
for file in *.sam
do
  sample=${file%%.sam}
  samtools view -Sb ${file} > ${sample}.bam
  rm -rf ${file}
done
```

### statistics

```
for file in *.bam
do
  sample=${file%%.bam}
  echo processing ${sample}...
  samtools view -c -f 4 ${file} > ${sample}.reads-unmapped.count.txt
  samtools view -c -F 4 ${file} > ${sample}.reads-mapped.count.txt
  samtools sort -T ${sample}.sorted -o ${sample}.sorted.bam ${file}
  samtools index ${sample}.sorted.bam
  samtools idxstats ${sample}.sorted.bam > ${sample}.reads.by.contigs.txt
done
```

### mapped/unmapped reads

```
echo "sample;mapped;unmapped" > mapping_stats.txt
for file in *.reads-mapped.count.txt
do
  sample=${file%%.reads-mapped.count.txt}
  mapped=`less ${file}`
  unmapped=`less ${sample}.reads-unmapped.count.txt`
  echo "${sample};${mapped};${unmapped}" >> mapping_stats.txt
done
```

## GENECALLING - FragGeneScan

### Link creation
`ln -s /home/ubuntu/miniconda3/pkgs/fraggenescan-1.31-hec16e2b_4/bin/train/ ./`

### run
`FragGeneScan -s final.contigs.fa -w 1 -o final.contigs_fgs -t complete -p 8`


## MAKE ANNOTATION TABLE

### make raw table

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Metagenome/Course_MB140P94/count-up-mapped-from-results-txt-with-ctg-length.py`

`python2.7 count-up-mapped-from-results-txt-with-ctg-length.py *.reads.by.contigs.txt`

### 2a NORMALISE MAPPING TABLE PER BASE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Metagenome/Course_MB140P94/normalize-mapping-table-by-read-length-and-ctg-length.py`

`python2.7 normalize-mapping-table-by-read-length-and-ctg-length.py summary-count-mapped.tsv 250 table_normalised.txt`

### 2b NORMALISE MAPPING TABLE PER SAMPLE

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Metagenome/Course_MB140P94/normalize_table_by_columns.py`

`python2.7 normalize_table_by_columns.py table_normalised.txt 2 1000000 table_normalised_per_sample.txt`

### MULTIPLY BY GENECALL

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Metagenome/Course_MB140P94/contig_mapping_to_genecall_mapping.py`

`python2.7 contig_mapping_to_genecall_mapping.py final.contigs_fgs.faa table_normalised_per_sample.txt`


## CONTINUE WITH ANNOTATION...

### get annotated metagenome

`wget -O SeqMe_assembly.zip "http://www.biomed.cas.cz/mbu/lbwrf/example/SeqMe_assembly.zip"`


`cd SeqMe_assembly/`

`wget https://raw.githubusercontent.com/pdobbler/cool-python-scripts/main/Metagenome/Course_MB140P94/link_simple_table_to_mapping_table.py`


`python2.7 link_simple_table_to_mapping_table.py MG_mapping_normalised_per_sample_genecall.txt TAXONOMY_BEST_OF_SIMPLE.txt TAX bitscore MG_NORM_TAX.tab`

`python2.7 link_simple_table_to_mapping_table.py MG_NORM_TAX.tab FOAM_KO_SIMPLE_MULTI.txt KEGG e-val MG_NORM_TAX_FOAM.tab`

`python2.7 link_simple_table_to_mapping_table.py MG_NORM_TAX_FOAM.tab CAZy_BEST_SIMPLE.txt CAZy e-val MG_NORM_TAX_FOAM_CAZy.tab`

```
mkdir FINAL_TABLE
mv MG_NORM_TAX_FOAM_CAZy.tab FINAL_TABLE/
mv *_tree.tab FINAL_TABLE/
cd FINAL_TABLE/
```


## ANNOTATION EXPLORER

`wget -O TABLE_EXPLORER_LINUX.zip "http://www.biomed.cas.cz/mbu/lbwrf/example/TABLE_EXPLORER_LINUX.zip"`

```
unzip TABLE_EXPLORER_LINUX.zip
cd TABLE_EXPLORER_LINUX/
bash run.sh
```

## binning contigs to genomes (BONUS)

```
for file in *.cut
do
  echo "${file}"
done > reads_list_file.txt
```

### run MaxBin2

`mkdir maxbin2_output`

`run_MaxBin.pl -contig final.contigs.fa -reads_list reads_list_file.txt -out maxbin2_output/maxbin2`




