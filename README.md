In this repository, I have added 2 projects of my BioInformatics course. Here is the brief Expplanation of each:

-------------------------------------------

# Read Mapping and Genome Assembly

## Overview

This project involves analyzing short-read sequencing data, performing quality control, creating genome assemblies, and conducting read mapping. The objective is to gain hands-on experience with bioinformatics tools and techniques for processing and analyzing whole-genome sequencing (WGS) data of *E. coli*.

## Project Structure

The project is divided into three main parts:

- **Part A:** Downloading E. coli WGS data, preliminary analyses, and quality controls.
- **Part B:** De novo genome assembly.
- **Part C:** Read mapping.

## Project Objectives

1. **Download and Analyze Short Reads:**
   - Download *E. coli* WGS data (SRR8185316) from the SRA database.
   - Perform preliminary analyses on the fastq files, including counting reads, extracting sequences, and plotting quality metrics.

2. **Quality Control:**
   - Perform quality control on the reads using FastQC.
   - Interpret and document the quality control results.

3. **De Novo Genome Assembly:**
   - Assemble the genome using SPAdes.
   - Evaluate the quality of the draft genome assembly using Quast and compare it to the reference genome.

4. **Read Mapping:**
   - Map the short reads to the reference genome using BWA.
   - Convert SAM files to indexed BAM files and visualize the mapped reads.
   - Analyze the read mapping results, including calculating the percentage of mapped reads and understanding the CIGAR format.


## Requirements

- **Tools:** SRA Toolkit, FastQC, SPAdes, Quast, BWA, Samtools, IGV
- **Languages:** Python for preliminary analyses and plotting.

--------------------------------------------------------------------------

# RNA-Seq Analysis

## Overview

This project involves RNA-Seq analysis of paired next-generation sequencing (NGS) expression profiles from normal and tumor colorectal tissues.
The study is based on the dataset with accession number GSE104836. The analysis includes quality control, read trimming,
read mapping, gene expression quantification, differential gene expression analysis, and Gene Ontology (GO) enrichment analysis.

## Project Structure

The project is divided into five main parts:

- **Part A:** Quality control and trimming of RNA-Seq data.
- **Part B:** Read mapping to the reference genome.
- **Part C:** Building a gene expression matrix.
- **Part D:** Differential gene expression analysis.
- **Part E:** Gene Ontology enrichment analysis.

## Data Preparation

We require to preprocess paired data of one patient. Convert the SRA files into forward and reverse Fastq.gz files using the SRA Toolkit.

## Part A: Quality Control and Trimming

1. **Quality Control:**
   - Assess the read qualities using FastQC.
   - Review the FastQC reports for any issues in the read quality.

2. **Trimming:**
   - Improve the read qualities using Trimmomatic.
   - Recheck the trimmed reads with FastQC to ensure quality improvement.

## Part B: Read Mapping

1. **Indexing the Reference Genome:**
   - Index the Homo sapiens (GRCh38) reference genome using HISAT2.
     

2. **Mapping Reads:**
   - Map reads to the reference genome using HISAT2.

    
3. **Converting SAM to BAM:**
   - Convert SAM files to BAM files and sort them using Samtools.


## Part C: Building Gene Expression Matrix

1. **Counting Aligned Reads:**
   - Use HTSeq to count aligned reads for differential expression analysis.

2. **Merging Count Files:**
   - Merge the individual count files into a single expression matrix.

## Part D: Differential Gene Expression Analysis

1. **Using edgeR:**
   - Perform differential gene expression analysis using the edgeR package in R.


## Part E: Gene Ontology Enrichment Analysis

1. **GOseq Analysis:**
   - Perform Gene Ontology enrichment analysis using the GOseq package in R.


## Requirements

- **Tools:** SRA Toolkit, FastQC, Trimmomatic, HISAT2, Samtools, HTSeq, edgeR, GOseq
- **Languages:** Python, Bash, and R for analysis and plotting
