#!/bin/bash

# Directories
KALLISTO_DIR=""
OUTPUT_FILE=""

# Initialize temporary file
TEMP_FILE="temp_counts.tsv"

# Find all SRR directories and sort them
SRR_DIRS=($(find "$KALLISTO_DIR" -maxdepth 1 -type d -name "SRR*" | sort))

# Loop through SRR directories, find abundance files, and extract counts
for SRR_DIR in "${SRR_DIRS[@]}"; do
    # Get the sample ID
    SRR_ID=$(basename "$SRR_DIR")
    # Find the abundance file
    ABUNDANCE_FILE=$(find "$SRR_DIR" -maxdepth 1 -type f -name "*_abundance.tsv" | head -n 1)
    # If this is the first file, extract the target_id column for the header
    if [[ -z "$FIRST_FILE" ]]; then
        cut -f1 "$ABUNDANCE_FILE" > "$TEMP_FILE"
        FIRST_FILE="$ABUNDANCE_FILE"
    fi

    # Extract the counts column and add to the temporary file
    cut -f4 "$ABUNDANCE_FILE" | tail -n +2 > "${SRR_ID}_counts.tmp"
    paste "$TEMP_FILE" "${SRR_ID}_counts.tmp" > "${TEMP_FILE}_new"
    mv "${TEMP_FILE}_new" "$TEMP_FILE"
    rm "${SRR_ID}_counts.tmp"
done

# Add header to the table
HEADER="Transcript_ID"
for SRR_DIR in "${SRR_DIRS[@]}"; do
    HEADER+="\t$(basename "$SRR_DIR")"
done
echo -e "$HEADER" > "$OUTPUT_FILE"
cat "$TEMP_FILE" >> "$OUTPUT_FILE"

# Remove temporary file
rm "$TEMP_FILE"
