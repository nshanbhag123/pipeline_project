#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 18:47:50 2023

@author: nirushanbhag
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:25:45 2023

@author: nirushanbhag
"""
import os
from os import popen
from Bio import Entrez
from Bio import SeqIO
import glob
import re

def SRR_donor(srr_no):
    my_dict = {"SRR5660030":"Donor 1 (2dpi)", "SRR5660033": "Donor 1 (6dpi)", "SRR5660044": "Donor 3 (2dpi)", "SRR5660045":"Donor 3 (6dpi)"}
    
    return my_dict[srr_no]
  
def build_index(main): 
    os.chdir(main)
    os.makedirs(main + '/bowtie_index')
    os.chdir(main + '/bowtie_index')

    
    handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype="fasta") 
    
    records = list(SeqIO.parse(handle, "fasta"))
    
    SeqIO.write(records, os.getcwd() + "/index.fasta", "fasta")
    
    bowtie_cmd1 = "bowtie2-build " + main + "/bowtie_index/index.fasta " + "HCMV"
    
    os.system(bowtie_cmd1)
    

def bowtie(main, input_reads_path):
    
    #changing the directory to the pipeline project folder
    os.chdir(main)
    
    #changing the directory to the trimmed fastq folder
    os.chdir(input_reads_path)
    
    #using glob to append the path names of files ending with fastq to a list to iterate through
    list_of_fastqs = []
    for name in glob.glob(input_reads_path + "/*.fastq"):
        list_of_fastqs.append(name)
        
    #sorting the list of fastqs
    list_of_fastqs = sorted(list_of_fastqs)
    
    #making a directory for the bowtie2 output
    os.makedirs(main + "/" + "bowtie2_output")
    
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
        bowtie_cmd2 = "bowtie2 " + "-x " + index + " -1 " + fastq1 + " -2 " + fastq2 + " -S " + main + "/bowtie2_output/" + SRR_no+ "_map.sam " + "--al-conc " + main + "/bowtie2_output/" + SRR_no + "_mapped_%.fq" 
        os.system(bowtie_cmd2)
        
        
        #changing directories into the output of bowtie
        os.chdir(main + "/bowtie2_output")
        

        #using grep to find the number of reads in the file after bowtie.There should be the same number of reads in each of the paired end mapped files, so I'm only looking at one of the output mapped files to get the number of reads. 
        
        end_reads = os.popen("grep -c '@SRR' " + SRR_no + "_mapped_1.fq").read().strip('\n')


        #writing the before and after reads to the log file
        log_file.write(SRR_donor(SRR_no) + " had " + start_reads + " read pairs before Bowtie2 filtering and " + end_reads + " read pairs after" + "\n")
        
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
        
        
        
        
def spades(main):
    os.chdir(main)
    os.makedirs(main + "/spades")
    os.chdir(main + "/spades")
    
    btlist = []
    for name in glob.glob(main + "/bowtie2_output/*.fq"):
        btlist.append(name)
        
    btlist = sorted(btlist)
    
    spades_cmd = "spades.py " + "-k 77, 99, 127 -t 8 --only-assembler " + "--pe-1 1 " + btlist[0] + " --pe-2 1 " + btlist[1] + " --pe-1 2 " + btlist[2] + " --pe-2 2 " + btlist[3] + " --pe-1 3 " + btlist[4] + " --pe-2 3 " + btlist[5] + " --pe-1 4 " + btlist[6] + " --pe-2 4 " + btlist[7] + " -o " + main + "/spades/" + "assembly/"
    
    log_file.write(spades_cmd + "\n")
      
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
        
    
    os.system(spades_cmd)
    
    
def calc_contigs(main):
    os.chdir(main + "/spades/assembly/")
    
    records = SeqIO.parse(main + "/spades/assembly/contigs.fasta", format = "fasta")
    
    no_contigs = 0
    contig_length = 0 
    longest_contig = ""
    for record in records:
        if len(record.seq) > len(longest_contig):
            longest_contig = record
        if len(record.seq) > 1000:
            no_contigs += 1
            contig_length += len(record.seq)
            
    
    log_file.write("There are " + str(no_contigs) + " contigs > 1000 bp in the assembly." + "\n")
    
    log_file.write("There are " + str(contig_length) + " bp in the assembly." + "\n")
    
    
          
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
        
    
    return(longest_contig)

def blast(main):
    os.makedirs(main + "/blast")
    os.chdir(main + "/blast")
    
    HCMV_db = open(main + "/blast/HCMV.fasta", "w")
    
    handle = Entrez.esearch(db = "nucleotide", term = ("Betaherpesvirinae[Organism] OR Betaherpesvirinae[All Fields]) AND refseq[filter]"), retmax = 50) #using esearch to find all the gen bank ids corresponding with betaherpesvirinae
    
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
    
    input_file = main + "/blast/HCMV.fasta"
    output_file = main + "/blast/HCMV_db"
    
    makeblastdb = "makeblastdb -in " + input_file + " -dbtype nucl " + "-out " + output_file
    
    os.system(makeblastdb)
    
    longest_contig = calc_contigs(main)
    
    query = open(main + "/blast/query.fasta", "w")
    
    query.write(">" + str(longest_contig.id))
    query.write("\n")
    query.write(str(longest_contig.seq))
    query.write("\n")
    
    query.close()
    
    input_file = main + "/blast/query.fasta"
    output_file = main + "/blast/HCMV.txt"
    db = main + "/blast/HCMV_db"
    
    blastn = "blastn -query " + input_file + " -db " + db + " -num_threads 4 -max_hsps 1 -max_target_seqs 10 -out " + output_file + " -outfmt '7 stitle pident length qstart qend sstart send bitscore evalue'"
    
    
    os.system(blastn)
    
    output_content = open(output_file, "r")
    
    log_file.write(output_content.read())
    
      
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
    log_file.write("\n")
    
    
input_reads_path = os.getcwd() + "/trimmed_fastq_files"
    
main = os.getcwd() + "/PipelineProject_Niru_Shanbhag"

os.makedirs(main)
os.chdir(main)

log_file = open(main + "/PipelineProject.log", "w")


build_index(main)
bowtie(main, input_reads_path)
spades(main)
blast(main)

log_file.close()
