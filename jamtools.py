"""
Written by Jack Duryea (duryea@bcm.edu)
Waterland Labs
Baylor College of Medicine
Children's Nutritional Research Center

MIT License

Copyright (c) 2017 Jack Duryea

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""


import sys
"""
This module provides a number of useful functions for COMP 182, including
manipulating graphs, plotting data, and timing functions.
"""

import matplotlib.pyplot as plt
import pylab
import types
import time
import math
import copy
import numpy as np
from collections import defaultdict
## Plotting functions

def show():
    """
    Do not use this function unless you have trouble with figures.

    It may be necessary to call this function after drawing/plotting
    all figures.  If so, it should only be called once at the end.

    Arguments:
    None

    Returns:
    None
    """
    plt.show()



def plot_lines(data, title="Title", xlabel="X", ylabel="Y", labels=None, filename=None):
    """
    Plot a line graph with the provided data.

    Arguments: 
    data     -- a list of dictionaries, each of which will be plotted 
                as a line with the keys on the x axis and the values on
                the y axis.
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    labels   -- optional list of strings that will be used for a legend
                this list must correspond to the data list
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    ### Check that the data is a list
    if not isinstance(data, types.ListType):
        msg = "data must be a list, not {0}".format(type(data).__name__)
        raise TypeError(msg)

    ### Create a new figure
    fig = pylab.figure()

    ### Plot the data
    if labels:
        mylabels = labels[:]
        for i in range(len(data)-len(labels)):
            mylabels.append("")
        for d, l in zip(data, mylabels):
            _plot_dict_bar(d, l)
        # Add legend
        pylab.legend(loc='best')
        gca = pylab.gca()
        legend = gca.get_legend()
        pylab.setp(legend.get_texts(), fontsize='medium')
    else:
        for d in data:
            _plot_dict_bar(d)

    ### Set the lower y limit to 0 or the lowest number in the values
    mins = [min(l.values()) for l in data]
    ymin = min(0, min(mins))
    pylab.ylim(ymin=ymin)

    ### Label the plot
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    ### Draw grid lines
    pylab.grid(True)

    ### Show the plot
    fig.show()

    ### Save to file
    if filename:
        pylab.savefig(filename)

def _dict2lists(data):
    """
    Convert a dictionary into a list of keys and values, sorted by
    key.  

    Arguments:
    data -- dictionary

    Returns:
    A tuple of two lists: the first is the keys, the second is the values
    """
    xvals = data.keys()
    xvals.sort()
    yvals = []
    for x in xvals:
        yvals.append(data[x])
    return xvals, yvals

def _plot_dict_line(d, label=None):
    """
    Plot data in the dictionary d on the current plot as a line.

    Arguments:
    d     -- dictionary
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if label:
        pylab.plot(xvals, yvals, label=label)
    else:
        pylab.plot(xvals, yvals)

def _plot_dict_bar(d, xmin=None, label=None):
    """
    Plot data in the dictionary d on the current plot as bars. 

    Arguments:
    d     -- dictionary
    xmin  -- optional minimum value for x axis
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if xmin == None:
        xmin = min(xvals) - 1
    else:
        xmin = min(xmin, min(xvals) - 1)
    if label:
        pylab.bar(xvals, yvals, align='center', label=label)
        pylab.xlim([xmin, max(xvals)+1])
    else:
        pylab.bar(xvals, yvals, align='center')
        pylab.xlim([xmin, max(xvals)+1])

def _plot_dict_scatter(d):
    """
    Plot data in the dictionary d on the current plot as points. 

    Arguments:
    d     -- dictionary

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    pylab.scatter(xvals, yvals)
    



# import argparse

# parser = argparse.ArgumentParser(description='Process some integers.')
# parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                     help='an integer for the accumulator')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')

# args = parser.parse_args()
# print args.accumulate(args.integers)

help_string = """
Hi there! Welcome to Jamtools, a 100% Python tool for analyzing Sam files. 
Here are some commands you can use:
Usage: python jamtools.py command [file]
count_reads:				counts the number of reads in a .Sam file
count_non_unique_mapping:		counts the number of reads that did not map uniquely i.e. have a MAPQ score of 0
count_duplicates:			counts the number of duplicate reads (PCR and/or optical)
count_unmapped_reads:			counts the number of reads that did not map
get_RONUM:				returns the rate of non unique mapping
"""

field_keys = ["QNAME", "FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL"]


# Checks to see if a file is in SAM format
def check_sam(file):
	if len(file) > 4 and file[-4:] != ".sam":
		print "error, file not in SAM format, check suffix"
		return False
	else:
		return True

# Checks to see if a line is a header in the file
def is_header(line):
	# Header lines begin with @
	return line[0] == "@"

# Takes a string representation of a decimal number and converts it to an 
# 11 bit binary number as a string
def to_binary(dec_str):
	decimal_value = int(dec_str)
	binary_value = bin(decimal_value)
	# Take off 0b
	binary_value = binary_value[2:]
	for i in range(11-len(binary_value)):
		binary_value = "0" + binary_value
	return binary_value



# Count the total number of reads in the sam file, ignore headers
def count_reads(samfile):
	"""
	Counts the total number of reads in the file
	Input: samfile - a file in SAM format
	Output: the number of reads
	"""

	# Make sure file is a sam file
	if not check_sam(samfile):
		return
	
	#Open file and read
	file = open(samfile, "r")
	count = long(0)
	for line in file:
		# Make sure we don't count headers
		if not is_header(line):
			count+=1
	print "number of reads:", count
	return count

# Count the reads that did not map uniquely
def count_non_unique_mapping(samfile):
	"""
	Counts the number of reads that have a MAPQ score of 0, and thus
	map to 2 or more locations with equal probability
	"""

	# Make sure file is a sam file
	if not check_sam(samfile):
		return

	# Open file
	file = open(samfile, "r")
	count = long(0)

	for line in file:
		# 
		if not is_header(line):
			read_data = {}
			# Split up the line into its fields
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# 0 indicates the read can map to 2 or more locations with equal probability
			if float(read_data["MAPQ"]) == 0:
				count += 1
	print count
	return count
	

# Count the number of PCR or optical duplicate reads
def count_duplicates(samfile):
	"""
	Counts the number of duplicate reads in the file. 
	Looks at the FLAG field of each line in the file, breaks
	it into its binary equivalent, and checks whether or not
	the bit 10 (first bit = bit 0) has been set. A 1 in bit 10 indicates that 
	a read is a duplicate from PCR or Illumina sequencing.
	"""
	
	# Make sure file is a sam file
	if not check_sam(samfile):
		return

	# Open file
	file = open(samfile, "r")
	count = long(0)

	for line in file:
		if not is_header(line):
			read_data = {}
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# Convert FLAG field to binary
			binary_value = to_binary(read_data["FLAG"])

			# if bit 10 is 1, then the read is a PCR duplicate
			if binary_value[0]=="1":
				count += 1
	print count
	return count

# Count the number of reads that did not map anywhere
def count_unmapped_reads(samfile):
	"""
	Count the number of reads that did not map, this
	is indicated in the SAM file under the FLAG field, 
	if bit 9 (first bit = bit 0) is set then the read did not map
	"""

	# Make sure file is a sam file
	if not check_sam(samfile):
		return
	
	file = open(samfile, "r")
	count = long(0)

	for line in file:
		if not is_header(line):
			read_data = {}
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# Convert FLAG field to binary
			binary_value = to_binary(read_data["FLAG"])

			# if bit 9 is 1, then the read is a PCR duplicate
			if binary_value[1] == "1":
				count += 1
	print count
	return count


# Returns the rate of non unique mapping for the samfile
# This is the number of reads that do not have MAPQ score of 42
# divided by the number of reads in total
# TODO: consider only reads that mapped
def get_RONUM(samfile):
	"""
	Counts the number of reads that have a MAPQ score of 0, and thus
	map to 2 or more locations with equal probability
	"""

	# Make sure file is a sam file
	if not check_sam(samfile):
		return

	# Open file
	file = open(samfile, "r")
	mapq_42_count = long(0)
	total_count = long(0)
	max_ronum = 0

	threshold = 42

	for line in file:
		# 
		if not is_header(line):
			total_count += 1
			read_data = {}
			# Split up the line into its fields
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# 0 indicates the read can map to 2 or more locations with equal probability
			if float(read_data["MAPQ"]) > max_ronum:
				max_ronum =  float(read_data["MAPQ"])
			if float(read_data["MAPQ"]) < threshold:
				mapq_42_count += 1

	print mapq_42_count/float(total_count)
	return mapq_42_count/float(total_count)


def get_mapping_efficiency_report(samfile):
	# Make sure file is a sam file
	if not check_sam(samfile):
		return

	# Open file
	file = open(samfile, "r")
	data = []
	rolling_sum = 0.0
	rolling_min = 100
	rolling_max = 0
	num_reads = long(0)
	sum_squared = long(0)
	data = defaultdict(lambda:0)
	for line in file:
		# 
		if not is_header(line):
			num_reads += 1

			read_data = {}

			# Split up the line into its fields
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# 0 indicates the read can map to 2 or more locations with equal probability
			score = int(read_data["MAPQ"])
			rolling_sum += score
			if score > rolling_max:
				rolling_max = score
			if score < rolling_min:
				rolling_min = score
			sum_squared += (score**2)
			data[score] += 1

	mean = rolling_sum/float(num_reads)
	plot_lines([data], filename = "MAPQ Scores")
	# Compute variance
	var = (sum_squared/num_reads) - (mean**2)
	sd = var**0.5
	print "max: ",rolling_max
	print "min: ",rolling_min
	print "mean: ",mean

	print "var: ",var
	print "sd: ", sd


# Reports a bunch of information about the SAM file
def full_report(samfile):
	return



# Creates a scatter plot of average RONUMS with varying sample sizes
# across the library
# TODO: make this process more random
# 
def library_complexity(samfile):
	# Make sure file is a sam file
	if not check_sam(samfile):
		return

	# Open file
	file = open(samfile, "r")
	num_reads = count_reads(samfile)

	# Partition into 100 data points
	num_points = 10
	bin_size = long(num_reads/num_points)
	data = {}
	non_unique_count = 0
	bin_num = 0
	bin_count = 0
	reads_processed = 0

	file = open(samfile, "r")
	for line in file:

		if bin_count >= bin_size:
			data[bin_num+1] = float(non_unique_count)/reads_processed
			bin_num += 1
			bin_count = 0

		# Make sure not a header
		if not is_header(line):
			reads_processed += 1
			bin_count += 1

			read_data = {}
			# Split up the line into its fields
			for key,value in zip(field_keys, line.split()):
				read_data[key] = value
			# 0 indicates the read can map to 2 or more locations with equal probability
			if float(read_data["MAPQ"]) == 0:
				non_unique_count += 1
	plot_lines([data], "RONUM", "Percentage of Library", "RONUM", labels=None, filename="RONUM_Plot")

# Entry, parse command line, use argparse later
if len(sys.argv) >= 2:
	command = sys.argv[1]
	if command == "help" or command == "h" or command == "--help":
		print help_string
	
	if command == "count_reads":
		if len(sys.argv) > 2:
			count_reads(sys.argv[2])

	if command == "count_non_unique_reads":
		if len(sys.argv) > 2:
			count_non_unique_reads(sys.argv[2])

	if command == "count_unmapped_reads":
		if len(sys.argv) > 2:
			count_unmapped_reads(sys.argv[2])

	if command == "get_RONUM":
		if len(sys.argv) > 2:
			get_RONUM(sys.argv[2])

	if command == "library_complexity":
		if len(sys.argv) > 2:
			library_complexity(sys.argv[2])

	if command == "get_mapping_efficiency_report":
		if len(sys.argv) > 2:
			get_mapping_efficiency_report(sys.argv[2])


