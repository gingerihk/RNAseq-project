#!/bin/bash

OUTPUT_DIR=""
QC_DIR=""

for fastq_file in "$OUTPUT_DIR"/*.fastq.gz; do
   fastqc -o "$QC_DIR" "$fastq_file"
done

#applied on trimmed file as well 

for trimmed_file in "$OUTPUT_DIR"/*.fastq.gz; do
   fastqc -o "$QC_DIR" "$trimmed_file"
done

