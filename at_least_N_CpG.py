#!/usr/bin/python
####################################################
# 
# Bins human genome into x-size bins and keeps track of the
# percentage of bins that have at least N CpG sites in them
#
# Jack Duryea
# Waterland lab, Baylor College of Medicine
# <duryea@bcm.edu>
####################################################


# Imports
import numpy as np
import sys
import pandas as pd
from os import listdir
from os.path import isfile, join
import os
import xlrd
from collections import defaultdict
import matplotlib.pyplot as plt
import pylab
import types
import collections
import glob

## Plotting function
def histogram(data, filename="AtLeastNCoverage.png", max_bins = 20):
    # data to plot
    n_groups = max_bins

    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.15
    opacity = 0.8

    
    rects1 = plt.bar(index - 2* bar_width, data[0][:max_bins], bar_width,
                     alpha=opacity,
                     color='b',
                     label='100bp')
    rects2 = plt.bar(index - 1* bar_width, data[1][:max_bins], bar_width,
                     alpha=opacity,
                     color='r',
                     label='200bp')
    rects3 = plt.bar(index, data[2][:max_bins], bar_width,
                     alpha=opacity,
                     color='g',
                     label='300bp')
    rects4 = plt.bar(index + 1* bar_width, data[3][:max_bins], bar_width,
                     alpha=opacity,
                     color='k',
                     label='400bp')
    rects5 = plt.bar(index +2* bar_width, data[4][:max_bins], bar_width,
                     alpha=opacity,
                     color='c',
                     label='500bp')
    

    hfont = {'fontname':'Arial'}
    plt.xlabel('Num CpGs sites',**hfont)
    plt.ylabel('Percentage of Bins with at least N cpg sites',**hfont)
    plt.xticks(index + bar_width, range(1, max_bins+1), **hfont)
    plt.tight_layout()
    
    # Higher resolution
    pylab.savefig(filename,format='png', dpi=1000)
    fig = plt.figure()
    plt.close(fig)
    #plt.show()
bin_size = 100 # Default bin size, vary this for experiments
genome_bins = {} # contains bin data
total_cpgs_in_genome = 0 # count total cpgs across genome
cpgs_per_chrome = {}



# Read each chromosome file
#files = os.listdir(os.getcwd())
files = ["Positions/chrX.CpG.positions.txt"]
files =  glob.glob("/Users/student/Desktop/MouseBinning/Positions/*")

for filename in files:
    if filename[-3:] == "txt": # If this is a text file
        chrome_name = filename[:filename.index(".")]
        my_file = open(filename,"r")
        bins = defaultdict(lambda:[])

        # Keeps track of cpgs
        all_cpgs_in_chrm = []
        for line in my_file:
            data = line.split() # line is tab seperated, first entry is chrm number, second is position
            cpg_pos = int(data[1])
            all_cpgs_in_chrm.append(cpg_pos)
            bin_num = int(cpg_pos/(bin_size+0.0)) # Divide by bin size to find bin number
            bins[bin_num].append(chrome_name+"-"+data[1])
            total_cpgs_in_genome += 1 # Count cpgs

        genome_bins[chrome_name] = bins
        cpgs_per_chrome[chrome_name] = len(all_cpgs_in_chrm) # Distribution of CpGs
    
## Now that we have the data, look at the coverage
num_cpgs_per_bin = defaultdict(lambda:0)
mouse_genome_size = 2730871774 # Number of DNA elements in the mm10 genome
total_number_of_bins = (mouse_genome_size)/bin_size
bins_per_size = defaultdict(lambda:0)

# Go through each chromosome
for chrm in genome_bins:
    for chrm_bin in genome_bins[chrm]: # Find bins
        num_cpgs = len(genome_bins[chrm][chrm_bin])
        bins_per_size[num_cpgs] += 1

data = []
for x in bins_per_size:
    data.append(bins_per_size[x])

plot_data = {}
for i in range(len(data)):
    plot_data[i+1] = (100.0*sum(data[i:]))/total_number_of_bins

histogram(plot_data, "at_leastN.png")
