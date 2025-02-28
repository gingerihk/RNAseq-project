# Co-expression and co-splicing analysis on raw RNA-Seq data 

## Introduction
This project aims to identify gene and transcript level co-expression and co-splicing patterns by utilizing WGCNA (Weighted Gene Co-Expression Network Analysis) R package. It is designed to be flexible and can be applied to any RNA-Seq dataset.

The project folder contains bash scripts used for downloading and processing raw data to generate a counts table and one R script that performs the desired analyses based on the counts table. 

## Overview

This workflow is designed to perform:
  - Data aquisition: downloading raw RNA-Seq data based on Project ID number
  - Pre-processing: QC and trimming
  - Pseudoalignement and transcript quantification
  - Co-expression and co-splicing analysis
  - Functional enrichment 

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
     
  2. **SRA-Toolkit** 
     - install SRA-Tools from git in your Software directory:
     `git clone https://github.com/ncbi/sra-tools.git\`
    
  3. **Conda**
     - HPCs should normally have this module installed (e.g. check it by running `which conda`);
     
  4. **FastQC, Cutadapt, Kallisto**
     - if installed on cluster, load them before running your command
     - if not installed on cluster, install them globally in your base environment or in your project environment.
     e.g. `conda install -c bioconda cutadapt`
  5. **Required R libraries**
     - install them by running `install.packages()` or `BiocManager::install()`
 
## Usage 
### job-submission.sh: submit and execute bash scripts on cluster 
  
### download-raw-rnaseq-data.sh: downloads FASTA files based on project ID

### fastqc.sh: performs quality check on each FASTQ file

### fastq-trimming.sh: trims 3' and 5' ends of each file
  
### kallisto.sh: explains how to build a transcriptome index and performs pseudoalignment 

### counts.sh: combines counts from multiple samples into a table
  
### RNA-Seq-project.R: processes the transcript counts table and performs co-expression and co-splicing analyses 

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
* Ensure all dependencies are loaded before running the scripts
* Adapt paths and module versions to your HPC environment 












  

  
