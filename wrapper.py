import os
import subprocess
from Bio import SeqIO
import pandas as pd

os.system('mkdir PipelineProject_Richa_Patel')
os.chdir("PipelineProject_Richa_Patel") 

os.system('gh repo clone richapatel138/PipelineProject_Richa_Patel')

log = open("PipelineProject.log", "w")

#PART 2

os.system('esearch -db nucleotide -query "NC_006273.2" | efetch -format gb > HCMV.gb')

with open("HCMV_CDS.fasta", "w") as fasta_out:
    record = SeqIO.read('HCMV.gb', "genbank")
    cds_count = 0

    for feature in record.features:
        if feature.type == "CDS":
            if "protein_id" in feature.qualifiers and "translation" in feature.qualifiers:
                protein_id = feature.qualifiers["protein_id"][0]
                sequence = feature.location.extract(record).seq
                
                fasta_out.write(f">{protein_id}\n{sequence}\n")
                cds_count += 1

log.write(f"The HCMV genome (NC_006273.2) has {cds_count} CDS.")
log.write("\n")

os.system('kallisto index -i index.idx HCMV_CDS.fasta')

#Part 3

os.system('mkdir results')

index_file = "index.idx"
output_dir = "results"
input_dir = "Seqs"
threads = 2
bootstraps = 30

with open('SRR_info') as f: #open the file
    srr_info = f.read().splitlines()

SRR = srr_info[0::3]
donor = srr_info[1::3]
conditions = srr_info[2::3]

#make sure output dict exits
os.makedirs(output_dir, exist_ok=True)

#Loop through each SRR id and run kallisto
for srr in SRR:
    output_path = os.path.join(output_dir, srr)
    input_file_1 = os.path.join(input_dir, f"{srr}_1.fastq")
    input_file_2 = os.path.join(input_dir, f"{srr}_2.fastq")

    print(f"Running kallisto for {srr}...")

    cmd = [
        "kallisto", "quant",
        "-i", index_file,
        "-o", output_path,
        "-b", str(bootstraps),
        "-t", str(threads),
        input_file_1,
        input_file_2
   ]

    os.system(cmd)

    print(f"Completed {srr}")

print("Done with Kallisto!")

#get TMP output
log.write(f'sample\tcondition\tmin_tpm\tmed_tpm\tmeam_tpm\tmax_tpm')
log.write("\n")
for item, name in zip(SRR, conditions):
    tsv_path = f'results/{item}/abundance.tsv'
    df = pd.read_csv(tsv_path, sep='\t')
    df_mean = df['tpm'].mean()
    df_min = df['tpm'].min()
    df_med = df['tpm'].median()
    df_max = df['tpm'].max()
    log.write(f'{item}\t{name}\t{df_min}\t{df_med}\t{df_mean}\t{df_max}\n')

#Part 4 - Sleuth
os.system('Rscript sleuth_commands.R')

import pandas as pd

df = pd.read_csv('sleuth_results.txt', sep=' ')
log.write("\n")
log.write(f"{df[['target_id', 'test_stat', 'pval', 'qval']]}")
log.write("\n")

#Part 5
os.system('datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report')
os.system('unzip ncbi_dataset.zip')
os.system('bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCNV')

index_file = "HCNV"

for srr in SRR:
    input_file_1 = f"{srr}_1.fastq"
    input_file_2 = f"{srr}_2.fastq"
    sam_file = f"{srr}map.sam"
    mapped_reads = f"{srr}_mapped_%.fq"

    print(f"Running Bowtie for {srr}...")

    cmd = [
        "bowtie2", 
        "-x", index_file,
        "-1", input_file_1,
        "-2", input_file_2,
        "-S", sam_file,
        "--al-conc", mapped_reads
   ]

    os.system(cmd)

    print(f"Completed {srr}")

print("Done with Bowtie!")

#number of reads before and after
#loop through each SRR id and calculate reads
for srr, donor, cond in zip(SRR,donor,conditions):
    input_file = f"{srr}_1.fastq"
    mapped_file = f"{srr}_mapped_1.fq"

    with open(input_file) as f: #open the file
        in_file = f.read().splitlines()

    with open(mapped_file) as f: #open the file
        map_file = f.read().splitlines()
    
    input = int(len(in_file) / 4)
    mapped = int(len(map_file) / 4)

    log.write(f'{donor} ({cond}) had {input} read pairs before Bowtie2 filtering and {mapped} read pairs after.\n')

#Part 6 - Spades Assembly

os.system("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq -o Donor_1_assembly/")
os.system("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660044_mapped_1.fq --pe-2 1 SRR5660044_mapped_2.fq --pe-1 2 SRR5660045_mapped_1.fq --pe-2 2 SRR5660045_mapped_2.fq -o Donor_3_assembly/")   

log.write("SPAdes Commands:\n")
log.write("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq -o Donor_1_assembly/")
log.write("\n")
log.write("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660044_mapped_1.fq --pe-2 1 SRR5660044_mapped_2.fq --pe-1 2 SRR5660045_mapped_1.fq --pe-2 2 SRR5660045_mapped_2.fq -o Donor_3_assembly/")
log.write("\n")

#Part 7 - BLAST
os.system('mkdir BLAST')
os.chdir("BLAST") 
os.system('datasets download virus genome taxon Betaherpesvirinae --refseq --include genome')
os.system('unzip ncbi_dataset.zip')
os.system('cd ..')
os.system('makeblastdb -in BLAST/ncbi_dataset/data/genomic.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl')

assembly_contigs = ['Donor_1_assembly/contigs.fasta', 'Donor_3_assembly/contigs.fasta']
first_contigs = ['Donor_1_contig', 'Donor_3_contig']

for assembly,contig in zip(assembly_contigs,first_contigs):
    with open(assembly) as fasta_file:
        first_seq = next(SeqIO.parse(fasta_file, "fasta"))
        with open(contig, "w") as output_file:
            SeqIO.write(first_seq, output_file, "fasta")

contigs_out_blast = ['Donor_1_BLAST', 'Donor_3_BLAST']

for contigs_in,blast_out in zip(first_contigs,contigs_out_blast):
    query_seqfile = contigs_in
    output_file = blast_out
    blast_command = "blastn -query " + query_seqfile+" -db betaherpesvirinae -max_hsps 1 -out "+output_file+" -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_command)

donor_1_output = subprocess.getoutput("head -10 Donor_1_BLAST")
donor_3_output = subprocess.getoutput("head -10 Donor_3_BLAST")

log.write(f'Donor1:\n')
log.write(f'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevaluet\tstitle\n')
log.write(donor_1_output)
log.write("\n")

log.write(f'Donor3:\n')
log.write(f'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevaluet\tstitle\n')
log.write(donor_3_output)
log.write("\n")

log.close()
