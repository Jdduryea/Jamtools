# Jamtools - Software for working with.SAM files
Welcome to Jamtools! A Python tool for working with SAM files created with SAM files (files resulting from mapped sodium bisulfite converted next generation sequencing data). Jamtools performs similar functions to the popular software tool Samtools, but simplifies the usage of many functions and improves upon the speed due Jamtool's lightweight nature. 

# Version
Version: 1.0
License: MIT

# Developer
Jamtools was created by J.D. Duryea (duryea@bcm.edu). Please feel free to ask any questions, voice concerns, or suggest improvements!

# Installation
Installation is simple: just download the jamtools.py file and copy it to your working directory.

# Functions
Currently, Jamtools offers only a handful of functions, but I'm constantly adding more! These function have been useful to us at Waterland labs and hopefully they can help you too!

## Usage
All functions in Jamtools have the same commandline usage form:
python jamtools.py example_command example.sam

count_reads:				        
counts the number of reads in a .Sam file

count_non_unique_mapping:		
counts the number of reads that did not map uniquely i.e. have a MAPQ score of 0

count_duplicates:			      
counts the number of duplicate reads (PCR and/or optical)

count_unmapped_reads:			  
counts the number of reads that did not map

get_RONUM:				          
returns the rate of non unique mapping (the number of reads that did not map uniquely/ the number of reads                                 total) A note about unique reads, the term "unique" is rather vague, we only have for each read a MAPQ score                               in the range [0, 42], which is the -10log_10 probability of how likely an the read is not placed correctly.                               High MAPQ score = good. For this function we just return the fraction of reads that have a MAPQ score of less                             than 42.

get_mapping_report:         
returns statistics about the MAPQ scores for the reads in the library. Saves a plot (png) showing the                                     distribution of scores.

library_complexity:         
Beta, creates a scatter plot of the percentage of the library analyzed vs the average RONUM of that sample.

full_report:                
Beta, returns a whole bunch of statistics about the library mapping, most of them coming from the functions                               about
