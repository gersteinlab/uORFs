#! /usr/bin/env python

# Author: Russell Ault

# uorf_gtf_to_GERP.py converts a uORF GTF file to the GERP file used by
# uorfs_GERPrate.py. This GERP format is basically the VAT file format
# with only a few fields

# python uorf_gtf_to_GERP.py CODON <uORF-gtf-file>

import os
import sys

global Codon	# start codon for the uORFs

# This function is no longer used
#def read_fields(fields):
#	return fields[0], int(fields[3]), int(fields[4]), 1


def codon_base_positions(fields_lst, strand='+', start=True, bases_to_extract=3):
	base_count = 0
	last_base = 0
	
	# default values for + strand and start codon
	chromPos = 3
	comparePos = 4
	nfield = 0
	
	if strand == '+' and not start:
		chromPos = 4
		comparePos = 3
		nfield = len(fields_lst) - 1
	if strand == '-':
		if start:
			chromPos = 4
			comparePos = 3
			nfield = 0
		else:
			chromPos = 3
			comparePos = 4
			nfield = len(fields_lst) - 1
			
	positions = []
	positions.append(int(fields_lst[nfield][chromPos]))
	base_count += 1
	while base_count < bases_to_extract:
		if positions[last_base] == int(fields_lst[nfield][comparePos]):
			#print "Reached a point in bases_codon where last field wasn't big enough"
			if start:
				nfield += 1
			else:
				nfield -= 1
			positions.append(int(fields_lst[nfield][chromPos]))
			base_count += 1
			last_base += 1
		else:
			if (strand == '+' and start) or (strand == '-' and not start):
				positions.append(positions[last_base] + 1)
			else:
				positions.append(positions[last_base] - 1)
			base_count += 1
			last_base += 1
		"""print positions[-1] # DEBUGGING CODE
		print "strand:\t%s" % strand"""
	return positions


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print "Usage: python uorfgtfToGERP.py CODON <uORF-gtf-file>"
		print "Incorrect number of arguments supplied"
		print "Exiting..."
		sys.exit(1)
	
	Codon = sys.argv[1]
	uORF_GTF_filename = sys.argv[2]
	
	uORF_GTF_file = open(uORF_GTF_filename, 'r')
	uORF_GERP_file = open("~/uORFs/uORFs_%s.gerp" % (Codon), "w")
	
	# Header line
	uORF_GERP_file.write('#chr\tstart\tlength\tstrand\tuORF-ID\tstart/end\n')
	
	prevUORF_ID = None
	fields_lst = []
	line = uORF_GTF_file.readline()  
	while line:
		fields = line.split()
		uORF_ID = fields[11].rstrip(';').strip('"').rstrip('.UTRCDSstop_codon')
		if uORF_ID != prevUORF_ID and prevUORF_ID:				# first line of a new uORF descriptor
			chromosome = fields_lst[0][0]
			strand = fields_lst[0][6]
			start_positions = codon_base_positions(fields_lst, strand, True)
			end_positions = codon_base_positions(fields_lst, strand, False)
			
			# Output to GERP file
			for base_position in start_positions:
				uORF_GERP_file.write('\t'.join([chromosome, str(base_position), '1', strand, prevUORF_ID, 'start']) + '\n')
			for base_position in end_positions:
				uORF_GERP_file.write('\t'.join([chromosome, str(base_position), '1', strand, prevUORF_ID, 'stop']) + '\n')
		
			# Reset fields list
			fields_lst = []

		fields_lst.append(fields)
		prevUORF_ID = uORF_ID
		line = uORF_GTF_file.readline()
		
		if not line:				# Process last uORF
			chromosome = fields_lst[0][0]
			strand = fields_lst[0][6]
			start_positions = codon_base_positions(fields_lst, strand, True)
			end_positions = codon_base_positions(fields_lst, strand, False)
			
			# Output to GERP file
			for base_position in start_positions:
				uORF_GERP_file.write('\t'.join([chromosome, str(base_position), '1', strand, prevUORF_ID, 'start']) + '\n')
			for base_position in end_positions:
				uORF_GERP_file.write('\t'.join([chromosome, str(base_position), '1', strand, prevUORF_ID, 'stop']) + '\n')
			
			break
	
	# Close files
	uORF_GTF_file.close()
	uORF_GERP_file.close()