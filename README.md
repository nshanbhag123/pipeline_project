# Python Pipeline Project

This pipeline project was part of an assignment for Computational Biology (COMP 383) in Spring 2023.

## Description

Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. Cheng et al. 2017
(https://www.ncbi.nlm.nih.gov/pubmed/29158406) sequences the transcriptomes of HCMV post infection. This python pipeline uses bioinformatics tools (Bowtie2, SPAdes, and command-line BLAST) to find the most similar HCMV strain to the two patient donor HCMV transcriptomes. 

This code was written for and tested in Python 3.8.

## Getting Started

### Software

* [SRA toolkit](https://github.com/ncbi/sra-tools)
* [Fastq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fastq-dump.html)
* [Bowtie2](https://github.com/BenLangmead/bowtie2)
* [SPAdes](https://github.com/ablab/spades)
* [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

### Packages
* os
* argsparse
* BioPython
* glob
* re

### Set Up Environment

* Use the following command to clone this repo to your machine
```
git clone https://github.com/nshanbhag123/pipeline_project.git
```
* Next, cd into the pipeline_project directory
```
cd pipeline_project
```
* To ensure that all the files have copied properly to your machine, list the contents of the pipeline_project directory to ensure all files and directories are present
```
ls 
```
* If all files and directories are present, you may run one of the two following commands. Make sure you are in the main pipeline_project directory, and not in a subdirectory. 


### Running the pipeline
* There are two ways you can run pipeline: running on the sample set (trimmed fastq files) and the big data set (untrimmed fastq files). Both output the results to a folder named "PipelineProject_Niru_Shanbhag"
* Running the big data set
 ```
 python3 genome_assembly.py
 ```
* Running the sample data
 ```
 python3 genome_assembly.py --sample
 ```
 If you want to run both commands, clone the repository in different directories so that the results do not get overwritten.
### Output Folders
* To view the contents of pipeline output, make sure to cd into the PipelineProject_Niru_Shanbhag
```
cd PipelineProject_Niru_Shanbhag
```
* PipelineProject.log
  * Contains a log file with a summary of the bowtie2 output, the complete spades command, and the blast output
* bowtie2_index
  * Contains the HCMV index files
* bowtie2_output
  * Contains the Bowtie2 mapped reads to the HCMV index
* spades
  * Contains the SPAdes assembly from all four transcriptomes
* blast
  * Contains the blast output from blasting the longest contig from the SPAdes assembly against the Betaherpesvirinae subfamily

