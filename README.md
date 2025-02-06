# RNA-Seq Analysis Pipeline

## Overview
This script automates RNA-Seq data analysis, from quality control to read alignment and quantification. It integrates several bioinformatics tools to streamline the process.

## Features
- **Quality Control & Trimming**: Uses `fastp` or `fastqc` + `trimmomatic`.
- **Read Alignment**: Uses `HISAT2` for aligning reads to a reference genome.
- **Quantification**: Uses `featureCounts` to generate a count matrix.
- **MultiQC Report**: Summarizes QC results in an interactive HTML report.

## Dependencies
Ensure the following software is installed before running the script:

- `fastp`
- `HISAT2`
- `samtools`
- `trimmomatic`
- `fastqc`
- `featureCounts`
- `multiqc`
- `Anaconda/Miniconda`

Installation:
```bash
conda install -c bioconda fastp hisat2 samtools trimmomatic fastqc subread multiqc
```

## Usage
### 1. Prepare Input Data
- Ensure raw FASTQ files are in the working directory.
- Files should follow the naming convention: `_R1_001.fastq.gz` and `_R2_001.fastq.gz` for paired-end reads.

### 2. Run the Script
```bash
bash RNASeqpipeline_x.sh
```
Follow on-screen prompts to:
- Set up dependencies
- Choose the number of processing threads
- Select an organism
- Provide or download genome index and annotation (GTF) files

### 3. Output Files
- **Processed FASTQ files**: Stored in `processed_files/`
- **Trimmed and QC reports**: Stored in `fastp_results/` or `fastqc_results/`
- **Aligned BAM files**: Stored in the working directory
- **FeatureCounts output**: Stored in `quants/all-featurecounts.txt`
- **MultiQC Report**: Generated in the working directory (`multiqc_report.html`)

## License
This script is licensed under **CC BY-NC 4.0**.

## Author
**Jaidev** (2025)

