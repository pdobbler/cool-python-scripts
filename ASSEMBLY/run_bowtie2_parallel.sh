#!/bin/bash

# Path to the Bowtie2 index
BOWTIE2_INDEX="Margaux_MG_Megahit_build/Margaux_MG_Megahit.build"

# Create output directory for BAM files
OUTPUT_DIR="bowtie2_results"
mkdir -p "$OUTPUT_DIR"

# Find unique file prefixes based on paired-end files
find . -name "*.pe.qc.fq.1.gz" | sed 's/.pe.qc.fq.1.gz//' > file_prefixes.txt

# Function to process each sample
process_sample() {
    PREFIX=$1
    OUTPUT_FILE="$OUTPUT_DIR/${PREFIX##*/}.sorted.bam"
    bowtie2 -x "$BOWTIE2_INDEX" \
        -1 "${PREFIX}.pe.qc.fq.1.gz" \
        -2 "${PREFIX}.pe.qc.fq.2.gz" \
        -U "${PREFIX}.se.qc.fq.gz" | \
    samtools view -bS - | \
    samtools sort -o "$OUTPUT_FILE" && \
    samtools index "$OUTPUT_FILE"
}

# Export function for parallel to use
export -f process_sample
export BOWTIE2_INDEX OUTPUT_DIR

# Run in parallel
parallel -j $(nproc) process_sample :::: file_prefixes.txt

# Cleanup
rm file_prefixes.txt

echo "All jobs completed. Results are in the '$OUTPUT_DIR' directory."


