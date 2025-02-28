# Co-expression and co-splicing analysis on raw RNA-Seq data 

## Introduction
This project aims to identify gene and transcript level co-expression and co-splicing patterns by utilizing WGCNA (Weighted Gene Co-Expression Network Analysis) R package. It is designed to be flexible and can be applied to any RNA-Seq dataset.

The project folder contains bash scripts used for downloading and processing raw data to generate a counts table and a R script that performs the desired analyses based on the counts table. 

## Overview

This workflow is designed to perform:
  - data aquisition: downloading raw RNA-Seq data based on Project ID number
  - pre-processing: QC and trimming
  - pseudoalignement and transcript quantification
  - co-expression and co-splicing analysis
  - functional enrichment 

## Prerequisites 

- High-performance computing (HPC) system
- R and RStudio 
- Basic knowledge of RNA-Seq data processing and R programming
- Understanding of co-expression network analysis (e.g. WGCNA)
- Conda: required to manage environments 

  ### Dependencies
- Edirect: required for accessing and querying NCBI's databases
- SRA-Toolkit: required for downloading RNA-Seq data from NCBI SRA
- FastQC: required for reads quality control
- Cutadapt: required for read processing stages
- Kallisto: required for pseudoalignment and transcript quantification
- R libraries (e.g. WGCNA, dplyr, ggplot) 


## Installation

  1. **EDirect Tools** 
     - install by running the following command in your Software directory:
     `sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)`
     - activate EDirect in your terminal session:
     `export PATH=$\{HOME\}/edirect:$\{PATH\}` (export in home)
  2. **SRA-Toolkit** 
     - install SRA-Tools from git in your Software directory:
     `git clone https://github.com/ncbi/sra-tools.git\`
     - check README.file to find download page for pre-built binaries and run the following command:
    `wget https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)`
     - unzip the tar.gz file you just downloaded
     `tar -xvzf sratoolkit.3.1.1-centos_linux64.tar.gz`
     - check folder for fasterq-dump command
     - run `fasterq-dump` in Bin to check if fasterq-dump is working
     - load sra-tools module before running fasterq-dump
     ** Optional ** - add fasterq-dump to your PATH environment to run it from any directory (no more need for module load)
         - go back to home directory and run `ls -a` to reveal hidden files
         - enter .bashrc using a command-line text editor (e.g. nano, vim) and paste the path to Bin 
         - run command source .bashrc to save changes; now fasterq-dump command should work no matter the directory
  3. **Conda**
     - HPCs should normally have this module installed (e.g. check it by running `which conda`);
     - if installed, load them by running the following command:
       `module load conda`
     - conda needs to have miniforge module loaded first 
  4. **FastQC, Cutadapt, Kallisto**
     - if installed on cluster, load them before running your command
     - if not installed on cluster, install them globally in your base environment or in your project environment.
         **For example**, for installing cutadapt in your base environment you would run:
         `conda install -c bioconda cutadapt`
         If you want to install it in your project environment, first load conda, then create your project environment by running the following command:
         `conda create -n rna-seq-env` (replace rna-seq-env with the name of your project)
         Then install cutadapt:
         `conda install -c bioconda cutadapt`
         Verify installation:
         `cutadapt --version`
         Your project environemnt can contain as many command-line programs as you need. Creating such an environment is useful when your project requires specific programs and versions. Moreover, it enables exact reproduction of your analysis by using the same tools versions every time.
         Activate/deactivate environment: `conda activate rna-seq-env` `conda deactivate`

  5. **Required R libraries**
     - install them by running `install.packages()` or `BiocManager::install()`
 
## Usage 
### run-jobs: job script sent to the cluster to run the bash scripts 
### download-raw-rnaseq-data.sh : downloads FASTA files based on project ID
`SRR_LIST=$(esearch -db sra -query "$PROJECT_ID" | efetch -format runinfo | cut -d',' -f1 | grep SRR)`
 * esearch searches sra and retrieves a list of unique identifiers that match the query
 * efetch retrived metadata table in csv format
 * cut and grep extract SRR numbers from the result
`FILTERED_SRR_LIST=$(echo "$SRR_LIST" | awk '$1 >= "SRR10260429" && $1 <= "SRR10260508"'| sort)
 echo "$FILTERED_SRR_LIST" > "$SCRATCH_DIR/CD14_SRR_list.txt"`
 * awk and sort filteres and sorts the SRR numbers in the given range
 * echo saves filtered list to a .txt file
`for SRR in $(cat "$SCRATCH_DIR/CD14_SRR_list.txt"); do
    echo "Downloading $SRR..."
    fasterq-dump --temp "$TEMP_DIR" --outdir "$OUTPUT_DIR" --split-files --details "$SRR"
    gzip "$OUTPUT_DIR/${SRR}"_*.fastq
done`
* cat reads the list
* fasterq-dump downloads FASTQ files and gzip compresses them; --split-files is chosen because of paired data

# run-fastqc.sh : performs quality check on each FASTQ file
`for fastq_file in "$OUTPUT_DIR"/*.fastq.gz; do
   fastqc -o "$QC_DIR" "$fastq_file"
done`

# fastq-trimming.sh : trims 3' and 5' ends of each file
`cutadapt -u 10 -u -20 -U 10 -U -20 -o "$TRIMMED_R1" -p "$TRIMMED_R2" "$R1" "$R2"`
* trims bases from both forward and reverse reads
# kallisto.sh : explains how to build a transcriptome index and performs pseudoalignment 
`kallisto quant -i "$INDEX_FILE" -o "$OUTPUT_DIR" -b 100 "$R1" "$R2"
if [[ -f "$OUTPUT_DIR/abundance.tsv" ]]; then
        mv "$OUTPUT_DIR/abundance.tsv" "$OUTPUT_DIR/${SRR}_abundance.tsv"
    fi`
* quant estimates transcript abundances through pseudoalignment
* -b means boostrap replicates will be performed
* renames abundance.tsv tables

# counts.sh : combines counts from multiple samples into a table
`SRR_DIRS=($(find "$KALLISTO_DIR" -maxdepth 1 -type d -name "SRR*" | sort))`
* finds all SRR directories within the Kallisto output directory
* sorts them alphabetically
  
`ABUNDANCE_FILE=$(find "$SRR_DIR" -maxdepth 1 -type f -name "*_abundance.tsv" | head -n 1)`
* finds abundance file in each SRR directory
  
`cut -f4 "$ABUNDANCE_FILE" | tail -n +2 > "${SRR_ID}_counts.tmp"
paste "$TEMP_FILE" "${SRR_ID}_counts.tmp" > "${TEMP_FILE}_new"
mv "${TEMP_FILE}_new" "$TEMP_FILE"
rm "${SRR_ID}_counts.tmp`
* cut extracts transcript counts from the 4th column for each sample
* paste combines the extracted counts with the existing temporary table
* mv updates the temporary counts table

`HEADER="Transcript_ID"
for SRR_DIR in "${SRR_DIRS[@]}"; do
    HEADER+="\t$(basename "$SRR_DIR")"
done
echo -e "$HEADER" > "$OUTPUT_FILE"
cat "$TEMP_FILE" >> "$OUTPUT_FILE"
`
* loops through the folder where SRR directories are, extracts their names and > creates the output file and writes the header row to it
* cat displays the content of the file that stores transcript counts and >> apends the count data to the output file

`rsync -xatv --bwlimit=5000  tosevsa2@hoffman2.idre.ucla.edu:/u/scratch/t/tosevsa2/table_counts.tsv /Users/andreeaiuhaniak/Desktop/software-project`
* rsync transfers files between remote and local systems; from source to destination

# RNA-Seq-project.R : processes the transcript_counts table and performs co-expression and co-splicing analyses 
* hashtags provide explanations for what lines do

## Input files:
* Raw FASTQ files
* Homo_sapiens.GRCh38.cdna.all.fa.gz from Ensembl to build reference transcriptome
* Homo_sapiens.GRCh38.110.gtf.gz GTF file from Ensembl
* GSE138746_series_matrix.txt metadata table
* f.1740434152.516.txt splicing factors file from SpliceAidF 

## Output files:
* Quality reports: fastQC HTML reports
* Quantification results: table_counts.tsv
* Co-expression Modules: WGCNA module assignments
* Figures: heatmaps, dendrograms, bar plots 

## Notes
* ensure all dependencies are loaded before running the scripts
* adapt paths and module versions to your HPC environment 












  

  
