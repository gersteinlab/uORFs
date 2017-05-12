#! /usr/bin/env python

# Author: Russell Ault

# make_positive_and_negative_set.py takes a uORF table and a list of positive data set uORF IDs and 
# saves two uORF table files, one for the positive uORFs and one for the unlabeled uORFs

# python make_negativeset.py  -d <uORF-data-table>[-p <positive-IDs-filename>]+

import os
import sys
from random import randint

from make_negativeset import *

if __name__ == '__main__':
	if len(sys.argv) < 2:
		ERROR("Usage: python make_negativeset.py  -d <uORF-data-table>[-p <positive-IDs-filename>]+")

	data_filename = None
	positive_filename = None
	data_lines = {}
	positive_IDs = set()
	neutral_IDs = set()
	keep_line = False
	
	i = 1
	while i < len(sys.argv):
		if sys.argv[i] == '-d':
			if (i+1) == len(sys.argv):
				ERROR("No data table provided")
			elif positive_filename:
				ERROR("Data table must come before positive table!")
			else:
				data_filename = sys.argv[i+1]
				data_lines = getTable(data_filename)
		elif sys.argv[i] == '-p':
			if (i+1) == len(sys.argv):
				ERROR("No positive uORF IDs provided")
			else:
				positive_filename = sys.argv[i+1]
				positive_IDs = getIDs(positive_filename, positive_IDs, positive=True)
		elif sys.argv[i] == '-r':
		 	neutral_filename = sys.argv[i+1]
		 	neutral_IDs = getIDs(neutral_filename, neutral_IDs, positive=True)
		else:
			print "Option provided not defined: %s" % (sys.argv[i])
			ERROR("Usage: python make_negativeset.py  -d <uORF-data-table>[-p <positive-IDs-filename>]+")
			
		i += 2
	if not data_filename or not positive_filename:
		ERROR("Data or positive set file name not provided")
	
	data_file = open(data_filename)
	positive_file = open(data_filename.split('.')[0] + '.positive', 'w')
	unlabeled_file = open(data_filename.split('.')[0] + '.unlabeled', 'w')
	neutral_file = open(data_filename.split('.')[0] + '.neutral', 'w')
	
	comment_line = data_file.readline()
	positive_file.write(comment_line)
	unlabeled_file.write(comment_line)
	neutral_file.write(comment_line)
	ID_pos = 0
	
	for line in data_file:
		if not line.startswith('#'):
			cleanedLine = line.strip()
			if cleanedLine:
				if line.split()[ID_pos] in positive_IDs:
					positive_file.write(line)
				elif line.split()[ID_pos] in neutral_IDs:
					neutral_file.write(line)
				else:
					unlabeled_file.write(line)
	data_file.close()
	positive_file.close()
	unlabeled_file.close()
	neutral_file.close()