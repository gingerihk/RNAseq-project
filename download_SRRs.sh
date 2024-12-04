!/bin/bash

# Define your SRR list file
SRR_LIST="$HOME/srr_list.txt"  # Path to the file with SRR IDs

# Get the SRR ID for this specific job using the task ID
SRR=$(sed -n "$\{SGE_TASK_ID\}p" $HOME/srr_list.txt)

# Define the output directory
OUTPUT_DIR="$SCRATCH/sra_data/fastq_files
TEMP_DIR="$SCRATCH/TEMP
echo "Downloading $SRR..."
for SRR in $SRR_LIST; do
    fasterq-dump --outdir "$OUTPUT_DIR" -t "$TEMP_DIR" -split--files --details $SRR\
    gzip $SRR*
done

#!/bin/bash

# Define your SRR list file
SRR_LIST="$HOME/srr_list.txt"  # Path to the file with SRR IDs

# Define the output directory
OUTPUT_DIR="$SCRATCH/sra_data/fastq_files"
TEMP_DIR="$SCRATCH/TEMP"

# Loop through each SRR ID in the list
while read -r SRR; do

    # Run fasterq-dump
    fasterq-dump --outdir "$OUTPUT_DIR" -t "$TEMP_DIR" --split-files --details "$SRR"

    # Compress the output FASTQ files
    gzip "$OUTPUT_DIR/${SRR}"_*

done < "$SRR_LIST"

#try 

