#!/bin/bash

# Directories
INPUT_DIR="/u/scratch/t/tosevsa2/sra_data/fastq_files"
OUTPUT_DIR="/u/scratch/t/tosevsa2/sra_data/trimmed_files"

# Process all paired-end files in the input directory
for R1 in "$INPUT_DIR"/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    # Extract base name for output files
    BASENAME=$(basename "$R1" _1.fastq.gz)

    # Define output file paths
    TRIMMED_R1="$OUTPUT_DIR/${BASENAME}_trimmed_1.fastq"
    TRIMMED_R2="$OUTPUT_DIR/${BASENAME}_trimmed_2.fastq"

    # Run Cutadapt for trimming
   ~/.conda/envs/cutadapt_env/bin/cutadapt -u 10 -u -20 -U 10 -U -20 -o "$TRIMMED_R1" -p "$TRIMMED_R2" "$R1" "$R2"

   gzip "$TRIMMED_R1"
   gzip "$TRIMMED_R2"
done   
