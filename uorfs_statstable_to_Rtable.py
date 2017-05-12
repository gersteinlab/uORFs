#! /usr/bin/env python

# Author: Russell Ault

# uorfs_statstable_to_Rtable.py converts my uorf statistics table to a format amenable
# to analysis in R. Mainly N_A fields are changed as well as nonnumerical fields for use in 
# machine learning algorithms.

# python uorfs_statstable_to_Rtable.py <stats_table>

# Output goes to standard out

import os
import sys

if len(sys.argv) != 2:
	print "Usage: python uorfs_statstable_to_Rtable.py <stats_table>"
	sys.exit(1)

# I should make the fields chosen a command line argument
last_field = -6
table = open(sys.argv[1])
print '\t'.join(table.readline().split('\t'))
for line in table:
	processed_line = line.replace('none', '0').replace('None', '0').replace('weak', '1').replace('strong', '2').\
			   			  replace('N_A', '0.0').split('\t')[:last_field]
	line_extension = line.split('\t')[last_field:]
	total_line = processed_line + line_extension
	print '\t'.join(total_line)
table.close()