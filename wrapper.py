import os
import subprocess
from Bio import SeqIO
import pandas as pd

os.system("mkdir PipelineProject_Richa_Patel") #make new directory that will have all output created
os.chdir("PipelineProject_Richa_Patel") #move into new directory

log = open("PipelineProject.log", "w") #create log file

#See steps for how to download the sequences and get paired-end reads on the Github under Part 1

''' PART 2: Get HCMV Genome
- Retrieve the HCMV genome (NC_006273.2) from NCBI.
- Extract coding sequences (CDS) into a FASTA file.
- Outputs CDS statistic to log file.
- Builds a kallisto index.
'''

#use esearch to get the nuceleotide genbank entry for accession number NC_006273.2, this will inclde all data for this entry
#store all of that data into a file named HCMV.gb
os.system('esearch -db nucleotide -query "NC_006273.2" | efetch -format gb > HCMV.gb')

#create a new file that has only the CDS info -- HCMV_CDS.fasta
#take the HCMV.gb file and search through it for the CDSs, and then determine the associated protein_id and sequence
#write that CDS id and sequence out to the HCMV_CDS.fasta file 
with open("HCMV_CDS.fasta", "w") as fasta_out:
    record = SeqIO.read('HCMV.gb', "genbank")
    cds_count = 0 #this keeps track for how many CDS sequences there are in the gb file, used to output CDS stat in oputput 

    for feature in record.features:
        if feature.type == "CDS":
            if "protein_id" in feature.qualifiers and "translation" in feature.qualifiers:
                protein_id = feature.qualifiers["protein_id"][0]
                sequence = feature.location.extract(record).seq
                
                fasta_out.write(f">{protein_id}\n{sequence}\n") #format fasta file so its protien id followed by seq
                cds_count += 1 #each time a sequence is written, add 1 to the CDS count

#write to the log file the number of CDSs we extracted 
log.write(f"The HCMV genome (NC_006273.2) has {cds_count} CDS.")
log.write("\n")
log.write("\n")

#create the kallisto index using the HCMV_CDS.fasta, which had the coding info. 
#The index file will be called 'index.idx'
os.system('kallisto index -i index.idx HCMV_CDS.fasta')


''' PART 3: Transcript Quantification with Kallisto
- Quantifies transcript abundance for provided sequence read files.
- Outputs TPM statistics to log file.
'''

#makes a dir that will store the kallisto results
os.system('mkdir results')

#specify kallisto flag input
index_file = "index.idx" #name of index file we created
output_dir = "results" #dir where we want the results
threads = 2 #number of threads
bootstraps = 30 #number of bootstraps 

#the SRR_info file lists the SRR (name of .fastq files), which donor the sample comes from (Donor 1 or 2), and the sample conditions (2dpi or 6dpi) for each transcriptome
with open('../SRR_info') as f: #open the SRR info file
    srr_info = f.read().splitlines() #store each line as a new item in the list

SRR = srr_info[0::3] #parse through list to seperate each type of into, each position in each list corresponds to the same transcriptome
donor = srr_info[1::3]
conditions = srr_info[2::3]
input_dir = "../" #specify where the SRR files are located

#make sure output dict exits
os.makedirs(output_dir, exist_ok=True)

#loop through each SRR id and run kallisto
for srr in SRR: 
    output_path = os.path.join(output_dir, srr) #make output path for specified srr 
    input_file_1 = os.path.join(input_dir, f"{srr}_1.fastq") #make input path for forward read for specified srr
    input_file_2 = os.path.join(input_dir, f"{srr}_2.fastq") #make input path for reverse read for specified srr

    print(f"Running kallisto for {srr}...")
    #specify the command with given input/output paths, run for each SRR
    cmd = [
        "kallisto", "quant",
        "-i", index_file,
        "-o", output_path,
        "-b", str(bootstraps),
        "-t", str(threads),
        input_file_1,
        input_file_2]

    subprocess.run(cmd, check=True)

    print(f"Completed {srr}")

print("done with Kallisto!")

#get TMP output
log.write(f'sample\tcondition\tmin_tpm\tmed_tpm\tmeam_tpm\tmax_tpm') #format the header in the .log
log.write("\n")
for item, name in zip(SRR, conditions): #for each SRR kallisto result
    tsv_path = f'results/{item}/abundance.tsv' #path to where the abundance.tsv is located
    df = pd.read_csv(tsv_path, sep='\t') #read in the file as a pd df
    df_mean = round(df['tpm'].mean(), 2) #calculate the mean, min, median, max
    df_min = round(df['tpm'].min(), 2)
    df_med = round(df['tpm'].median(), 2)
    df_max = round(df['tpm'].max(), 2)
    log.write(f'{item}\t{name}\t{df_min}\t{df_med}\t{df_mean}\t{df_max}\n') #write mean, min, median, max along with identifier info to .log


''' PART 4: Differential Expression Analysis with Sleuth
- Uses Sleuth (R package) for differential expression analysis.
- Outputs significant transcripts (FDR < 0.05) results to log file.
'''
#run the sleuth commands to create the slueth results
os.system('Rscript ../sleuth_commands.R')

#open the sleuth results as a pd df
df = pd.read_csv('sleuth_results.txt', sep=' ')
log.write("\n")
log.write(f"{df[['target_id', 'test_stat', 'pval', 'qval']]}") #write out the sleuth results to the .log file
log.write("\n")
log.write("\n")

'''PART 5: Genome Alignment with Bowtie2
- Downloads and indexes a reference genome for HCMV.
- Aligns sequencing reads to the reference genome.
- Logs the number of read pairs before and after filtering.
'''
#dowload the HCMV refrence genome NCBI accession NC_006273.2 using datasets
os.system('datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report')
os.system('unzip ncbi_dataset.zip') #unzip the file
os.system('bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCNV') #build a bowtie index file (refrence)

index_file = "HCNV" #name of index file for bowtie mapping
input_dir = "../" #directory where the SRR reads exist

#for each SRR
for srr in SRR:
    input_file_1 = os.path.join(input_dir, f"{srr}_1.fastq") #make input path for forward reads for specified srr
    input_file_2 = os.path.join(input_dir, f"{srr}_2.fastq") #make input path for reverse reads for specified srr
    sam_file = f"{srr}map.sam" #specify name of .sam file
    mapped_reads = f"{srr}_mapped_%.fq" #specify name of mapped reads

    print(f"Running Bowtie for {srr}...")

    #specify the command with given input/output names, run for each SRR
    #this will just keep the reads that map to the refrence genome 
    cmd = [
        "bowtie2", 
        "-x", index_file,
        "-1", input_file_1,
        "-2", input_file_2,
        "-S", sam_file,
        "--al-conc", mapped_reads]

    subprocess.run(cmd, check=True)

    print(f"Completed {srr}")

print("Done with Bowtie!")

#number of reads before and after
#loop through each SRR id and calculate reads before and after bowtie mapping
for srr, donor, cond in zip(SRR,donor,conditions): #use the lists that we name earlier
    input_file = os.path.join(input_dir, f"{srr}_1.fastq") #path to raw data
    mapped_file = f"{srr}_mapped_1.fq" #path to mapped data

    with open(input_file) as f: #open the file
        in_file = f.read().splitlines() #read raw data

    with open(mapped_file) as f: #open the file
        map_file = f.read().splitlines() #read mapped data
    
    input = int(len(in_file) / 4) #calculate raw read number (divide by 4 since each entry takes 4 lines)
    mapped = int(len(map_file) / 4) #calculate mapped read number (divide by 4 since each entry takes 4 lines)

    #write number of input/mapped reads to .log file in specified format
    log.write(f'{donor} ({cond}) had {input} read pairs before Bowtie2 filtering and {mapped} read pairs after.\n') 

'''PART 6: De Novo Assembly with SPAdes
- Runs SPAdes for de novo genome assembly of mapped reads.
- Outputs assembled contigs for each donor sample.
- Writes SPAdes command to log file.
'''

#spades.py --> running spades 
# -k 77 --> specify k-mer length, this was given to us in assignment details
# -t 2 --> threads to use 
# --only-assembler --> only assemble
# Since we have 4 files for each donor (forward/reverse at 2dpi and 6dpi), specify the paired end sets
# --pe-1 1 SampleSRR5660030_mapped_1.fq --> 2dpi forward reads donor 1
# --pe-2 1 SampleSRR5660030_mapped_2.fq --> 2dpi reverse reads donor 1
# --pe-1 2 SampleSRR5660033_mapped_1.fq --> 6dpi forward reads donor 1
# --pe-2 2 SampleSRR5660033_mapped_2.fq --> 6dpi reverse reads donor 1
# -o Donor_1_assembly/ --> specify directory name for spades output
#same format was used for donor 3

#run assembly for donor 1
os.system("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SampleSRR5660030_mapped_1.fq --pe-2 1 SampleSRR5660030_mapped_2.fq --pe-1 2 SampleSRR5660033_mapped_1.fq --pe-2 2 SampleSRR5660033_mapped_2.fq -o Donor_1_assembly/")
#run assembly for donor 3
os.system("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SampleSRR5660044_mapped_1.fq --pe-2 1 SampleSRR5660044_mapped_2.fq --pe-1 2 SampleSRR5660045_mapped_1.fq --pe-2 2 SampleSRR5660045_mapped_2.fq -o Donor_3_assembly/")   

#write the spades commands used to generate assemblies to .log 
log.write("\n")
log.write("SPAdes Commands:\n")
log.write("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SampleSRR5660030_mapped_2.fq --pe-1 2 SampleSRR5660033_mapped_1.fq --pe-2 2 ampleSRR5660033_mapped_2.fq -o Donor_1_assembly/")
log.write("\n")
log.write("spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660044_mapped_1.fq --pe-2 1 SampleSRR5660044_mapped_2.fq --pe-1 2 SampleSRR5660045_mapped_1.fq --pe-2 2 SampleSRR5660045_mapped_2.fq -o Donor_3_assembly/")
log.write("\n")
log.write("\n")

'''PART 7: BLAST Analysis
- Downloads all herpesvirus viral genomes from NCBI.
- Creates a BLAST database.
- Extracts and BLASTs the longest contig from each donor.
- Logs the top 10 BLAST hits for each donor.
'''
#name new directory where ref seqs for database curation will be downloaded
os.system('mkdir BLAST') #name new directory where ref seqs for database curation will be downloaded
os.chdir("BLAST") #move into BLAST dir 
os.system('datasets download virus genome taxon Betaherpesvirinae --refseq --include genome') #download all reseq Betaherpesvirinae seqs
os.system('unzip ncbi_dataset.zip') #unzip the file
os.chdir("..") #move back into previous dir
os.system("cp BLAST/ncbi_dataset/data/genomic.fna genomic.fna") #cp the .fna file from the BLAST dir (contains the seq used to build ref database) to working dir
os.system('makeblastdb -in genomic.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl') #create blast db using Betaherpesvirinae seqs, title it betaherpesvirinae

os.system("cp Donor_1_assembly/contigs.fasta 1contigs.fasta") #copy the contigs created from Donor_1_assembly to working directory 
os.system("cp Donor_3_assembly/contigs.fasta 3contigs.fasta") #copy the contigs created from Donor_3_assembly to working directory 

assembly_contigs = ['1contigs.fasta', '3contigs.fasta'] #specify input .fasta files to get first contigs (longest contig)
first_contigs = ['Donor_1_contig', 'Donor_3_contig'] #specify the output file where we will write the first contig

for assembly,contig in zip(assembly_contigs,first_contigs): #for each input/ouput pair
    with open(assembly) as fasta_file: #from the input file
        first_seq = next(SeqIO.parse(fasta_file, "fasta")) #extract the first contig
        with open(contig, "w") as output_file: #write to the oputput file 
            SeqIO.write(first_seq, output_file, "fasta") #the extracted contig in .fasta format

contigs_out_blast = ['Donor_1_BLAST', 'Donor_3_BLAST'] #specify name for BLAST results file

for contigs_in,blast_out in zip(first_contigs,contigs_out_blast): #for each contig/output pair
    query_seqfile = contigs_in #the query file is the first contig we extracted
    output_file = blast_out #oputput file name will be what we specified in list
    #format the blast command which seraches in the betaherpesvirinae db we created, only keep the best alignment for any single query-subject pair of sequences
    #blast results will include Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title.
    blast_command = "blastn -query " + query_seqfile+" -db betaherpesvirinae -max_hsps 1 -out "+output_file+" -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_command) #run blast command

#from the blast results for each donor, get only the top 10 and write to a file
donor_1_output = subprocess.getoutput("head -10 Donor_1_BLAST")
donor_3_output = subprocess.getoutput("head -10 Donor_3_BLAST")

#format output file with specified reqs 
log.write(f'Donor1:\n')
log.write(f'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevaluet\tstitle\n')
log.write(donor_1_output) #write the top 10 blast results from file
log.write("\n")
log.write("\n")

#do the same with Donor 3
log.write(f'Donor3:\n')
log.write(f'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevaluet\tstitle\n')
log.write(donor_3_output)
log.write("\n")

#close the log as we are done!
log.close()
