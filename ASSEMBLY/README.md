### MT+MG ASSEMBLY

### CAZy

```
## using DBCAN3 for cazy annotation::

# split fasta #
python2.7  /mnt/DATA1/priscila/scripts_metagenome/split_fasta_by_group_size.py  /mnt/DATA/projects/priscila/bartek_reads_for_ass/gene_prediction/all_predicted_metagenemark_renamed2_2.faa 83000

mkdir SPLIT 
mv *.fas SPLIT
cd SPLIT
### run dbcan3
conda activate run_dbcan3

for file in *.fas
do
 sample=${file%%.fas}
 mkdir ${sample}
done

for file in *.fas
do
  sample=${file%%.fas}
  echo "run_dbcan ${file} protein --db_dir /mnt/DATA_FastRun0/DATABASES/run_dbcan3_database/db/ -t hmmer --out_dir ${sample} --hmm_cpu 1"
done > dbcan.sh
cat dbcan.sh | parallel

echo "" > all_dbCAN.txt
for file in *.fas
do
 sample=${file%%.fas}
 wc -l ${sample}/hmmer.out
 cat ${sample}/hmmer.out >> all_dbCAN.txt
done



# dbCAN ANNOTATION  --- this keeps only one hit per 'gene' - some genes have multiple hits to either binding domains and/or catalytic domains!!


export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k3,3 -k5,5g all_dbCAN.txt | sort -u -k3,3 --merge > all_dbCAN_best.txt
awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt | sort | uniq > hmm_names_uniq.txt
awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt > hmm_names.txt
awk -F'\t' '{print $3"\t"$5}' all_dbCAN_best.txt > all_dbCAN_best_gene_eval.txt
paste -d"\t" all_dbCAN_best_gene_eval.txt hmm_names.txt > CAZy_BEST_SIMPLE.txt


####### USE THIS to get MULTI CAZy per gene
dbcan_CAZy_multi_simple.py
python dbcan_CAZy_multi_simple.py all_dbCAN.txt all_dbCAN_simple_MULTI.txt

## add to bit table

python2.7 /mnt/DATA1/priscila/scripts_metagenome/link_simple_table_to_mapping_table.py BY_GENE_NORMALISED_per_GENECALL_TAX_LCA_SSU_ITS_MAG_CAZy_KOfam_KOG.tab all_dbCAN_simple_MULTI.txt CAZyme e-value BY_GENE_NORMALISED_per_GENECALL_TAX_LCA_SSU_ITS_MAG_CAZy_KOfam_KOG_multiCAZY.tab
```



### MAPPING

`chmod +x run_bowtie2_parallel.sh`

`./run_bowtie2_parallel.sh`

### COUNT (UN)MAPPED READS

```
# Output file for the summary table
OUTPUT_TABLE="reads_summary_table.txt"

# Add header to the output table
echo -e "SAMPLE_NAME\tMAPPED_READS\tUNMAPPED_READS" > "$OUTPUT_TABLE"

# Loop through BAM files and calculate reads
for file in *.sorted.bam; do
    sample=${file%%.sorted.bam}
    
    # Count unmapped and mapped reads
    unmapped=$(samtools view -c -f 4 "$file")
    mapped=$(samtools view -c -F 4 "$file")

    # Append the results to the table
    echo -e "${sample}\t${mapped}\t${unmapped}" >> "$OUTPUT_TABLE"

    # Reads by contig file
    samtools idxstats ${file} > ${sample}.reads.by.contigs.txt
done

echo "Summary table created: $OUTPUT_TABLE"
```