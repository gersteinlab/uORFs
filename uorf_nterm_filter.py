#! /usr/bin/env python

# Author: Russell Ault

# python uorf_nterm_filter.py [-u | -n] <uORF_stats_table> <uORFs.fa>

# Options:
#	-u : prints out the uORFs, whether totally within the 5' UTR or overlapping the CDS
#	-n : prints out the N-terminal extensions

# Note: The same Codon stats table and FASTA uORF file must be used together

# Output goes to standard output

from make_negativeset import ERROR

import os
import sys


def make_dict(filename, start_string="chr", key_pos=2, value_pos=5):
	file = open(filename)
	dict = {}
	
	for line in file:
		if line.startswith(start_string):
			dict[line.split()[0].split('|')[key_pos]] = line.split()[0].split('|')[value_pos]
			
	file.close()
	return dict
	
def print_stats_line(line, uORFs_dict, flag, ID_pos=0):
	classif = uORFs_dict.get(line.split()[ID_pos])
	if not classif:
		ERROR("uORFs dictionary was made incorrectly")
	elif flag == '-u':
		if classif != "CDSfull":
			print line,
	elif flag == '-n':
		if classif == "CDSfull":
			print line,
	else:
		ERROR("Invalid Flag provided")

if __name__ == '__main__':
	if len(sys.argv) != 4:
		ERROR("Usage: python uorf_nterm_filter.py [-u | -n] <uORF_stats_table> <uORFs.fa>")
	
	i = 1
	if not sys.argv[i].startswith('-'):
		ERROR("First argument must be a -u | -n flag")
	flag = sys.argv[i]
	uORF_stats_filename = sys.argv[i+1]
	uORFs_fasta_filename = sys.argv[i+2]
	
	uORFs_dict = make_dict(uORFs_fasta_filename)
	uORF_stats_file = open(uORF_stats_filename)
	print uORF_stats_file.readline(),			# Prints header line
	for line in uORF_stats_file:
		print_stats_line(line, uORFs_dict, flag)
	
	uORF_stats_file.close()
