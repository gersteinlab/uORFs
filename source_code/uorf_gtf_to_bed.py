#! /usr/bin/env python

# Author: Russell Ault
# uorf_gtf_to_bed converts a uORF GTF file to a BED file line by line. Thus, each exonic portion
# of the uorf is on its own line.

# How to call
# uorf_gtf_to_bed.py CODON

# Output goes to file names of the form: uORF_CODON.bed

import os
import sys

global uORFbedFile
#global sys.stdin
global Codon

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print "Incorrect number of arguments; need a codon argument!"
		print "Exiting..."
		sys.exit(1)
		
	Codon = sys.argv[1]
	
	sys.stdin = open("~/uORFs/uORFs_" + Codon + ".gtf", "r")
	uORFbedFile = open("~/uORFs/uORF_" + Codon + ".bed", "w")
	
	uORFbedFile.write("#chrom\tchromStart\tchromEnd\tname\tfake-score\t" + \
	                  "strand\n")
	
	# Now read through each upstream open reading frame, converting the GTF file to BED format.
	line = sys.stdin.readline()
	while line:
	    fields = line.split()
	    uORF_ID = fields[11].rstrip(';').strip('"').rstrip('.UTRCDSstop_codon')
	    chromosome = fields[0]
	    chromStart = str(int(fields[3]) - 1)
	    chromEnd = fields[4]
	    strand = fields[6]
	    score = 0		# This is a place holder for the BED format to allow the strand
	    	
	    # Output to BED file
	    uORFbedFile.write(''.join(['\t'.join([chromosome, chromStart, chromEnd, uORF_ID, str(score), strand]),'\n']))
	    
	    line = sys.stdin.readline()
	
	sys.stdin.close()
		
	
