# PipelineProject (Richa Patel)

Dependencies: 
- Entrez NCBI (https://www.ncbi.nlm.nih.gov/books/NBK179288/) 
- Kallisto (https://pachterlab.github.io/kallisto/manual)
- NCBI Datasets (https://github.com/ncbi/datasets)
- Bowtie2 (https://github.com/BenLangmead/bowtie2)
- Spades (https://github.com/ablab/spades)
- BLAST+ (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- Biopython (https://biopython.org/)
- fasterq-dump (https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- Please see the links for installation and documentation.

Packages: 
- os (import os)
- subprocess (import subprocess)
- pandas (import pandas as pd)
- SeqIO (from Bio import SeqIO)

## Pipeline Steps

**NOTE: Part 1 Needs to be done locally on the machine**

### Part 1: Retrieve Transcriptome Sequences
Part 1: Getting the Seqs
- Download wanted sequences from NCBI
  - Command format: wget [link]
    - Donor 1 (2dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030`
    - Donor 1 (6dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033`
    - Donor 3 (2dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044`
    - Donor 3 (6dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045`
- Convert to paired end fastq files
  - Command format: fasterq-dump SRRXXX
    - Donor 1 (2dpi): `fasterq-dump SRR5660030`
    - Donor 1 (6dpi): `fasterq-dump SRR5660033`
    - Donor 3 (2dpi): `fasterq-dump SRR5660044`
    - Donor 3 (6dpi): `fasterq-dump SRR5660045`

**NOTE: Parts 2-7 are automated though the wrapper python script.**

### Part 2: Get HCMV Genome
- Retrieve the HCMV genome (NC_006273.2) from NCBI.
- Extract coding sequences (CDS) into a FASTA file.
- Outputs CDS statistic to log file.
- Builds a kallisto index.

### Part 3: Transcript Quantification with Kallisto
- Quantifies transcript abundance for provided sequence read files.
- Outputs TPM statistics to log file.

### Part 4: Differential Expression Analysis with Sleuth
- Uses Sleuth (R package) for differential expression analysis.
- Outputs significant transcripts (FDR < 0.05) results to log file.

### Part 5: Genome Alignment with Bowtie2
- Downloads and indexes a reference genome for HCMV.
- Aligns sequencing reads to the reference genome.
- Logs the number of read pairs before and after filtering.

### Part 6: De Novo Assembly with SPAdes
- Runs SPAdes for de novo genome assembly of mapped reads.
- Outputs assembled contigs for each donor sample.
- Writes SPAdes command to log file.

### Part 7: BLAST Analysis
- Downloads all herpesvirus viral genomes from NCBI.
- Creates a BLAST database.
- Extracts and BLASTs the longest contig from each donor.
- Logs the top 10 BLAST hits for each donor.

## Creating Necessary Input Files

### "SRR_info"
- Lists the SRR (name of .fastq files), which donor the sample comes from (Donor 1 or 2), and the sample conditions (2dpi or 6dpi) for each transcriptome
- Each piece of info exist on a new line, repeates in same order for each transcriptome

### "sleuth_table.txt"
- Formatted to contain the sample name, condition, and path to kallisto output
- Needs to be space or tab seperated, needs a header, needs to end with a newline

## Running with Test Input
- In order or run the wrapper script with test input, first clone this repository to whatever directory you are working in:
  - `git clone https://github.com/richapatel138/COMP483_PipelineProject.git`
- Move into the cloned respository
  - `cd COMP483_PipelineProject`
- Run the wrapper .py script
  - `python wrapper.py`
- There should be a new directory (`PipelineProject_Richa_Patel`) with all created output and the .log file.
  - `cd PipelineProject_Richa_Patel`
- To view the contents of the log file:
  - `cat PipelineProject.log`


