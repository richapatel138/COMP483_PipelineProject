# PipelineProject (Richa Patel)

Dependencies: 
- Entrez NCBI (https://www.ncbi.nlm.nih.gov/books/NBK179288/) 
- Kallisto (https://pachterlab.github.io/kallisto/manual)
- NCBI Datasets (https://github.com/ncbi/datasets)
- Bowtie2 (https://github.com/BenLangmead/bowtie2)
- Spades (https://github.com/ablab/spades)
- BLAST+ (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- Biopython (https://biopython.org/)
- fasterq-dump (https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump_

Packages: 
- os (import os)
- subprocess (import subprocess)
- pandas (import pandas as pd)

## Part 1 (Retrieve Transcriptomes) 
Part 1: Getting the Seqs
1. Download wanted sequences from NCBI
- Command format: wget [link]
  - Donor 1 (2dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030`
  - Donor 1 (6dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033`
  - Donor 3 (2dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044`
  - Donor 3 (6dpi): `wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045`
2. Convert to paired end fastq files
- Command format: fasterq-dump SRRXXX
  - Donor 1 (2dpi): `fasterq-dump SRR5660030`
  - Donor 1 (6dpi): `fasterq-dump SRR5660033`
  - Donor 3 (2dpi): `fasterq-dump SRR5660044`
  - Donor 3 (6dpi): `fasterq-dump SRR5660045`

## Part 2 (Quantify TPM) 
1. Build transcriptome index for HCMV (NCBI accession NC_006273.2)
2. Output the number 




