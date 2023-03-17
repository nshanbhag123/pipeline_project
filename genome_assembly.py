from os import popen
from Bio import Entrez
from Bio import SeqIO
import os
import argparse
import glob
import re

def convert_sra_to_fastq(main):
    os.chdir(main)
    
    #make a directory of untrimmmed fastq files and cd into it
    os.makedirs(main + "/untrimmed_fastqs") 
    os.chdir(main + "/untrimmed_fastqs")
    
    #fastq dump command
    fastq_dump = "/Users/nirushanbhag/Downloads/Software/sratoolkit.3.0.0-mac64/bin/fastq-dump -I --split-files "

    #finding path to sra files
    for sra in glob.glob(main + '/SRAfiles/*'):
        fastqdump_cmd = fastq_dump + sra 
        os.system(fastqdump_cmd)

def get_files(main):
    os.chdir(main)

    #make a directory of SRA files and cd into it
    os.makedirs(main + '/SRAfiles') #make a directory to the SRA files
    os.chdir(main + '/SRAfiles') #change directory to the sra file folder path
    
    #getting the SRA files using wget
    donor1 = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030"
    donor12 = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033"
    donor3 = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044"
    donor32 = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"
    
    donor_list = [donor1, donor12, donor3, donor32]
    
    for k in donor_list: #retreiving the sra files for each donor
        os.system("/usr/local/bin/wget " + k)

def SRR_donor(srr_no):
    #function that matches SRR number with Donor Name. Given number, will give you name
    my_dict = {"SRR5660030":"Donor 1 (2dpi)", "SRR5660033": "Donor 1 (6dpi)", "SRR5660044": "Donor 3 (2dpi)", "SRR5660045":"Donor 3 (6dpi)"}
    
    return my_dict[srr_no]
  
def build_index(main): 
    os.chdir(main)
    
    #making a directory of bowtie index files and cd into it
    os.makedirs(main + '/bowtie_index')
    os.chdir(main + '/bowtie_index')
    
    #Using efetch to create the HCMV fasta, parsing it, and creating the index fasta file
    handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype="fasta") 
    records = list(SeqIO.parse(handle, "fasta"))
    SeqIO.write(records, main + "/bowtie_index" + "/index.fasta", "fasta")
    
    #Using bowtie2-build to create an index out of the index fasta file
    bowtie_cmd1 = "/Users/nirushanbhag/Downloads/Software/bowtie2-2.2.0/bowtie2-build " + main + "/bowtie_index/index.fasta " + "HCMV"
    
    os.system(bowtie_cmd1)
    
def bowtie(main, input_reads_path):
    os.chdir(main)
    
    #changing the directory to the input reads path (could be trimmed or untrimmed)
    os.chdir(input_reads_path)
    
    #using glob to append the path names of files ending with fastq to a list
    list_of_fastqs = []
    for name in glob.glob(input_reads_path + "/*.fastq"):
        list_of_fastqs.append(name)
        
    #sorting the list of fastqs
    list_of_fastqs = sorted(list_of_fastqs)
    
    #making a directory for the bowtie2 output
    os.makedirs(main + "/" + "bowtie2_output")
    
    #writing the header for the Bowtie2 output
    log_file.write("# Bowtie2 output" + "\n")
    
    for i in range(0, len(list_of_fastqs), 2): #iterating through the fastq list in pairs
        
        fastq1 = list_of_fastqs[i] #first fastq 
        fastq2 = list_of_fastqs[i+1] #second fastq
        
        #extracting the SRR number from the fastq using reg expressions
        match = re.search("SRR\d+", fastq1) 
        SRR_no = str(match.group()) 
        
        #using grep to find the number of reads in the file before bowtie
        start_reads = os.popen("grep -c '@SRR' " + fastq1).read().strip('\n') #
        
        #path to the index
        index = main + "/bowtie_index/HCMV"
        
        #running bowtie on the paired end reads
        bowtie_cmd2 = "/Users/nirushanbhag/Downloads/Software/bowtie2-2.2.0/bowtie2  " + "-x " + index + " -1 " + fastq1 + " -2 " + fastq2 + " -S " + main + "/bowtie2_output/" + SRR_no+ "_map.sam " + "--al-conc " + main + "/bowtie2_output/" + SRR_no + "_mapped_%.fq" 
        os.system(bowtie_cmd2)
        
        #changing directories into the output of bowtie
        os.chdir(main + "/bowtie2_output")
    
        #using grep to find the number of reads in the file after bowtie. There should be the same number of reads in each of the paired end mapped files, so I'm only looking at one of the output mapped files to get the number of reads. 
        
        end_reads = os.popen("grep -c '@SRR' " + SRR_no + "_mapped_1.fq").read().strip('\n')

        #writing the before and after reads to the log file
        log_file.write(SRR_donor(SRR_no) + " had " + start_reads + " read pairs before Bowtie2 filtering and " + end_reads + " read pairs after" + "\n")

    log_file.write("\n" + "\n" + "\n" + "\n")
        
def spades(main):
    os.chdir(main)
    
    #making a new directory for spades output
    os.makedirs(main + "/spades")
    
    #using glob to iterate through all files ending with fastq and adding those paths to a list
    btlist = []
    for name in glob.glob("**/*.fq", recursive = True):
        btlist.append(name)

    #sort the list 
    btlist = sorted(btlist)
    
    #paired end spades assembly
    spades_cmd = "/Users/nirushanbhag/Downloads/Software/SPAdes-3.15.4-Darwin/bin/spades.py " + "-k 77, 99, 127 -t 8 --only-assembler " + "--pe-1 1 " + btlist[0] + " --pe-2 1 " + btlist[1] + " --pe-1 2 " + btlist[2] + " --pe-2 2 " + btlist[3] + " --pe-1 3 " + btlist[4] + " --pe-2 3 " + btlist[5] + " --pe-1 4 " + btlist[6] + " --pe-2 4 " + btlist[7] + " -o " + "spades/" + "assembly/"
    
    #writing the spades command to the log file
    log_file.write("# Spades Assembly " + "\n")
    log_file.write(spades_cmd)
    log_file.write("\n" + "\n" + "\n" + "\n")
        
    #executing the spades command
    os.system(spades_cmd)
    
def calc_contigs(main):
    os.chdir(main + "/spades/assembly/")
    
    #parsing the contigs fasta from spades assembly
    records = SeqIO.parse(main + "/spades/assembly/contigs.fasta", format = "fasta")
    
    #code to find the longest contig and how many contigs greater than 1000 bp were in the assembly
    no_contigs = 0
    contig_length = 0 
    longest_contig = ""
    for record in records:
        if len(record.seq) > len(longest_contig):
            longest_contig = record
        if len(record.seq) > 1000:
            no_contigs += 1
            contig_length += len(record.seq)
            
    #writing contig information to the log file
    log_file.write("# Summary stats from SPAdes assembly output" + "\n")
    log_file.write("There are " + str(no_contigs) + " contigs > 1000 bp in the assembly." + "\n")
    
    log_file.write("There are " + str(contig_length) + " bp in the assembly." + "\n")
    
    log_file.write("\n" + "\n" + "\n" + "\n")
        
    
    return(longest_contig) #return the longest contig

def blast(main):
    os.chdir(main)
    
    #making a blast directory and cd into it
    os.makedirs(main + "/blast")
    os.chdir(main + "/blast")
    
    #preparing to write to the HCMV db, which we will blast our longest contig against
    HCMV_db = open(main + "/blast/HCMV.fasta", "w")
    
    handle = Entrez.esearch(db = "nucleotide", term = ("Betaherpesvirinae[Organism] OR Betaherpesvirinae[All Fields])"), retmax = 20000) #using esearch to find all the gen bank ids corresponding with betaherpesvirinae
    
    record = Entrez.read(handle) #reading the handle
    
    list_of_ids = record["IdList"] #extracting the ids
    
    handle = Entrez.efetch(db = "nucleotide", id = list_of_ids, rettype = "fasta") #fetching the actual information about the ids from the ids we got from esearch
    
    records = list(SeqIO.parse(handle, format = "fasta")) #parsing the handle
    
    
    for record in records: #for each record, write it to the multifasta file for HCMV db
        HCMV_db.write(">" + str(record.description))
        HCMV_db.write("\n")
        HCMV_db.write(str(record.seq))
        HCMV_db.write("\n")
        
    HCMV_db.close() #close the file
    
    input_file = main + "/blast/HCMV.fasta" #input file is the fasta file containing all Betaherpesvirinae subfamily sequences 
    output_file = main + "/blast/HCMV_db" #output location and name
    
    #make blast db command
    makeblastdb = "/Users/nirushanbhag/Downloads/Software/ncbi-blast-2.13.0+/bin/makeblastdb -in " + input_file + " -dbtype nucl " + "-out " + output_file
    
    os.system(makeblastdb) #running the makeblastdb command
    
    longest_contig = calc_contigs(main) #running the calc_contigs function to find the longest contig
    
    #opening a file called query fasta and writing to it the longest contig's id and sequence.
    query = open(main + "/blast/query.fasta", "w")
    query.write(">" + str(longest_contig.id))
    query.write("\n")
    query.write(str(longest_contig.seq))
    query.write("\n")
    query.close()
    
    #for the blastn, the input file is the query fasta (longest contig from assembly)
    input_file = main + "/blast/query.fasta"
    output_file = main + "/blast/HCMV.txt"
    db = main + "/blast/HCMV_db"
    
    #blast n command, with outfmt 6
    blastn = "/Users/nirushanbhag/Downloads/Software/ncbi-blast-2.13.0+/bin/blastn -query " + input_file + " -db " + db + " -num_threads 4 -max_hsps 1 -max_target_seqs 10 -out " + output_file + " -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    
    #executing the blastn command
    os.system(blastn)
    
    #writing the blast header and blast output to the log file
    output_content = open(output_file, "r")
    log_file.write("# BLAST output" + "\n")
    log_file.write("sacc" + "\t" + "pident" + "\t" + "length" + "\t" +  "qstart" + "\t" + "qend" + "\t" + "sstart" + "\t" + "send" + "\t" + "bitscore" + "\t" + "evalue" + "\t" + "stitle" + "\n"
)
    log_file.write(output_content.read())
  
## RUN FUNCTIONS 

#initializing arg parse. You can either run the pipeline on the trimmed dataset or the untrimmed dataset
parser = argparse.ArgumentParser()
parser.add_argument("--sample", help="run the pipeline on the trimmed fastq dataset", action = "store_true")
args = parser.parse_args()

#using the cwd as the parent directory of the pipeline output folder
main = os.getcwd() + "/PipelineProject_Niru_Shanbhag"
os.makedirs(main)

#if the user specifies the sample dataset
if args.sample:
    input_reads_path = os.getcwd() + "/trimmed_fastq_files"#input reads are in trimmed fastq files
    
    log_file = open(main + "/PipelineProject.log", "w") #opening the log file

    #Running the functions
    build_index(main)
    bowtie(main, input_reads_path)
    spades(main)
    blast(main)
    
#else run the untrimmed dataset
else:
    input_reads_path = main + "/untrimmed_fastqs" #input reads are in the untrimmed fastq files
    log_file = open(main + "/PipelineProject.log", "w") #opening the log file

    #Running the  data retreival functions
    get_files(main)
    convert_sra_to_fastq(main)

    #Running the bowtie2, spades, and blast functions
    build_index(main)
    bowtie(main, input_reads_path)
    spades(main)
    blast(main)
    
log_file.close() #close log file
