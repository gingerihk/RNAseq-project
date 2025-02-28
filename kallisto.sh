 Build index
-first download fasta file: wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O /u/scratch/t/tosevsa2/kallisto/Homo_sapiens.GRCh38.cdna.all.fa.gz
-then unzipp it bc kallisto needs it unzipped: gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz 
-build index: kallisto index -i /u/scratch/t/tosevsa2/kallisto/transcriptome.idx /u/scratch/t/tosevsa2/kallisto/Homo_sapiens.GRCh38.cdna.all.fa


#!/bin/bash
# run multiple jobs: visualize and check each joblog file created

# Directories
INPUT_DIR=""
OUTPUT_DIR=""
INDEX_FILE=""

# Create a list of unique SRR IDs by processing only `_1.fastq.gz` files
SRR_LIST=($(ls "$INPUT_DIR" | grep "_trimmed_1.fastq.gz" | sed 's/_trimmed_1.fastq.gz//' | sort))

# Get the SRR ID corresponding to this job
SRR=${SRR_LIST[$((SGE_TASK_ID-1))]}  # Adjust for job array index (1-based)

# Define input files
R1="$INPUT_DIR/${SRR}_trimmed_1.fastq.gz"
R2="$INPUT_DIR/${SRR}_trimmed_2.fastq.gz"

# Define sample name and output directory

# Run Kallisto
~/.conda/envs/kallisto/bin/kallisto quant -i "$INDEX_FILE" -o "$OUTPUT_DIR" -b 100 "$R1" "$R2"

# Rename abundance.tsv to include the sample name
    if [[ -f "$OUTPUT_DIR/abundance.tsv" ]]; then
        mv "$OUTPUT_DIR/abundance.tsv" "$OUTPUT_DIR/${SRR}_abundance.tsv"
    fi
