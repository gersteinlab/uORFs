#! /usr/bin/env python

# Author: Russell Ault

# python uorf_gtf_to_PhyloCSF.py <uORF.gtf file>

# Output goes to standard output

# I am getting trouble with multiple field uORFs, cutting off the last field when I shouldn't

from uorf_gtf_to_GERP import codon_base_positions
from make_negativeset import ERROR

import os
import sys

global STARTPOS
global ENDPOS

STARTPOS = 3
ENDPOS = 4


def DEBUG_FIELDS(uORF_fields, msg='', nfields=False):
	if nfields and nfields > 1:
		print msg
		for fields in uORF_fields:
			print '\t'.join(fields)
		

def get_uORF_lines(GTF_file, ID_pos=11, strand_pos=6):
	prevUORF_ID = None
	prevstrand = None
	fields_lst = []
	line = GTF_file.readline()  
	while line:
		fields = line.split()
		uORF_ID = fields[ID_pos].rstrip(';').strip('"').rstrip('.UTRCDSstop_codon')
		strand = fields[strand_pos]
		if uORF_ID != prevUORF_ID and prevUORF_ID:				# first line of a new uORF descriptor
			yield fields_lst, prevUORF_ID, prevstrand
			fields_lst = []
			
		fields_lst.append(fields)
		prevUORF_ID = uORF_ID
		prevstrand = strand
		line = GTF_file.readline()
		
		if not line:						# Process last uORF
			yield fields_lst, prevUORF_ID, prevstrand
			break


def in_interval(a, b, n):
	if a <= n <= b or a >= n >= b:
		return True
	else:
		return False

def should_remove_last_field(uORF_fields, last_uORF_base):
	last_field = -1		# same for both strands
	if in_interval(int(uORF_fields[last_field][STARTPOS]), int(uORF_fields[last_field][ENDPOS]),
				   last_uORF_base):
		return False
	else:
		return True


def remove_last_field(uORF_fields):
	return uORF_fields[:-1]

	
def modify_last_field(uORF_fields, last_uORF_base, strand):
	new_fields = uORF_fields[:-1]
	if strand == '+':
		new_fields.append(uORF_fields[-1][:ENDPOS] + [str(last_uORF_base)] + 
											 uORF_fields[-1][ENDPOS+1:])
	else:
		new_fields.append(uORF_fields[-1][:STARTPOS] + [str(last_uORF_base)] + 
											 uORF_fields[-1][ENDPOS:])
	return new_fields 
		
			
def trim_stop_codon(uORF_fields, strand):
	end_positions = codon_base_positions(uORF_fields, strand, False, 4)			# this works fine, I thionk
	#print end_positions
	end_position = -1 				# same for both strands
	last_uORF_base = end_positions[end_position]
	while should_remove_last_field(uORF_fields, last_uORF_base):
		#print "This ID gave a remove last field:\t%s" % (uORF_fields[0][11])
		uORF_fields = remove_last_field(uORF_fields)
	uORF_fields = modify_last_field(uORF_fields, last_uORF_base, strand)
	if strand == '-':
		uORF_fields.reverse()
	return unite_blocks(uORF_fields, strand)


def unite_blocks(uORF_fields, strand):
	if len(uORF_fields) == 1:
		return uORF_fields
	new_fields = []
	blocks_united = 0
	i = 0
	while (i+1) < len(uORF_fields):
		if int(uORF_fields[i][ENDPOS]) == int(uORF_fields[i+1][STARTPOS]) - 1:
			new_fields.append(uORF_fields[i][:STARTPOS+1] + uORF_fields[i+1][ENDPOS:])
			blocks_united += 1
			i += 1
		else:
			new_fields.append(uORF_fields[i])
		i += 1
	if (i+1) == len(uORF_fields):
		new_fields.append(uORF_fields[i])		# append last field
	
	if blocks_united > 0:
		return unite_blocks(new_fields, strand)		# iterates over and over until all blocks have been united
	else:
		if strand == '-':
			new_fields.reverse()
		return new_fields
		
	#ERROR("I have to debug unite blocks!")		
		
		# I am going to need to iterate this until now change is done

def get_uORF_intervals(uORF_fields, chrom=0):
	intervals = sorted([field[STARTPOS] + '-' + field[ENDPOS] for field in uORF_fields], key=lambda x: x.split('-')[0])
	chromosome = uORF_fields[0][chrom]
	delim = '+' + chromosome + ':'
	outputString = delim.join(intervals)
	outputString = chromosome + ':' + outputString
	return outputString
	
def print_to_PhyloCSF(uORF_ID, strand, uORF_intervals):
	print '\t'.join([uORF_intervals, strand, uORF_ID])	
			
if __name__ == '__main__':
	if len(sys.argv) != 2:
		ERROR("Usage: python uorf_gtf_to_PhyloCSF.py <uORF.gtf file>")
	
	GTF_file = open(sys.argv[1], 'r')
	strand_pos = 6
	uORFs = get_uORF_lines(GTF_file)
	uORF_fields, uORF_ID, strand = next(uORFs, [False, False, False])
	"""print uORF_fields
	print uORF_ID
	print strand
	#ERROR("Debugging")"""
	while uORF_fields:
		"""new_uORF_fields = uORF_fields.reverse()
		print uORF_fields
		print new_uORF_fields"""
		#print "THIS IS THE START OF THE UORF_ID:\t%s" % uORF_ID
		#DEBUG_FIELDS(uORF_fields, "This is before doing anything", len(uORF_fields))
		if strand == '-':
			uORF_fields.reverse() 
		uORF_fields = unite_blocks(uORF_fields, strand)
		#DEBUG_FIELDS(uORF_fields, "This was after unite_blocks", len(uORF_fields))
		uORF_fields = trim_stop_codon(uORF_fields, strand)
		#DEBUG_FIELDS(uORF_fields, "This was after trimming the stop codon (where I think my error is)", len(uORF_fields))
		uORF_intervals = get_uORF_intervals(uORF_fields)
		#DEBUG_FIELDS(uORF_fields, "This was after getting the intervals", len(uORF_fields))
		print_to_PhyloCSF(uORF_ID, strand, uORF_intervals)
		uORF_fields, uORF_ID, strand = next(uORFs, [False, False, False])
	GTF_file.close()
