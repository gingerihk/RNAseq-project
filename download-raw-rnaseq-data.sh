#!/bin/bash

PROJECT_ID="PRJNA577035"
SCRATCH_DIR=""
OUTPUT_DIR=""
QC_DIR=""
TEMP_DIR="$SCRATCH/TEMP"

# Search SRR numbers associated with PRJNA577035
SRR_LIST=$(esearch -db sra -query "$PROJECT_ID" | efetch -format runinfo | cut -d',' -f1 | grep SRR)

# Save the full SRR list to a file
echo "$SRR_LIST" | sort > "$SCRATCH_DIR/SRR_list.txt"

# Filter CD14
FILTERED_SRR_LIST=$(echo "$SRR_LIST" | awk '$1 >= "SRR10260429" && $1 <= "SRR10260508"'| sort)

# Save the filtered SRR list to a file in scratch
echo "$FILTERED_SRR_LIST" > "$SCRATCH_DIR/CD14_SRR_list.txt"

# Download and zip files
for SRR in $(cat "$SCRATCH_DIR/CD14_SRR_list.txt"); do
    echo "Downloading $SRR..."
    fasterq-dump --temp "$TEMP_DIR" --outdir "$OUTPUT_DIR" --split-files --details "$SRR"
    gzip "$OUTPUT_DIR/${SRR}"_*.fastq
done
