#!/usr/bin/python
####################################################
# 
# Bins human genome into x-size bins and keeps track of the
# percentage of CpG sites covered
#
# Jack Duryea
# Waterland lab, Baylor College of Medicine
# <duryea@bcm.edu>
####################################################
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

## Plotting function
def histogram(data,filename):
    # data to plot
    n_groups = len(data)
    data1 = data.values()
    #data2 = non_cancer.values()

    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.8

    
    rects1 = plt.bar(index + bar_width, data1, bar_width,
                     alpha=opacity,
                     color='b',
                     label='CpGs')
    
    # rects2 = plt.bar(index + 1.5*bar_width, data1, bar_width,
    #                  alpha=opacity,
    #                  color='#2188CD',
    #                  label='Cancer')

    
    hfont = {'fontname':'Arial'}
    plt.xlabel('Num CpGs sites',**hfont)
    plt.ylabel('Percentage of CpGs covered genome wide',**hfont)
    plt.xticks(index + bar_width, range(1, len(data)+1), **hfont)

    plt.tight_layout()
    # Higher resolution
    pylab.savefig(filename,format='png', dpi=500)
    fig = plt.figure()
    plt.close(fig)
    #plt.show()


bin_size = 200 # Default bin size, vary this for experiments
genome_bins = {} # contains bin data
total_cpgs_in_genome = 0 # count total cpgs across genome
cpgs_per_chrome = {}

# Read each chromosome file
for filename in os.listdir(os.getcwd()):
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

# Go through each chromosome
for chrm in genome_bins:
    for chrm_bin in genome_bins[chrm]: # Find bins
        num_cpgs = len(genome_bins[chrm][chrm_bin])
        num_cpgs_per_bin[num_cpgs] += num_cpgs

data = []
for x in num_cpgs_per_bin:
    data.append(num_cpgs_per_bin[x])

plot_data = {}
for i in range(len(data)):
    plot_data[i+1] = (100.0*sum(data[i:]))/total_cpgs_in_genome

print plot_data

histogram(plot_data, "bin200Data.png")


