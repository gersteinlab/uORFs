#! /usr/bin/env python

# Author: Russell Ault

# python all_positive_dataset.py <literature_positive_data_set_file>

import sys
import os

if len(sys.argv) != 2:
	print "Usage: python all_positive_dataset.py <literature_positive_data_set_file>"
	print "Exiting..."
	sys.exit(1)

positive_data = open(sys.argv[1])

for line in positive_data:
	fields = line.split()
	for ID in fields:
		print ID

positive_data.close()