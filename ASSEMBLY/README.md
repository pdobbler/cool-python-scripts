### MT+MG ASSEMBLY

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

    # Reads by contig file
    samtools idxstats ${file} > ${sample}.reads.by.contigs.txt

    # Append the results to the table
    echo -e "${sample}\t${mapped}\t${unmapped}" >> "$OUTPUT_TABLE"
done

echo "Summary table created: $OUTPUT_TABLE"
```
