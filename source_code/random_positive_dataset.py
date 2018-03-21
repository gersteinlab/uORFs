#! /usr/bin/env python

# Author: Russell Ault

# python random_positive_dataset.py <literature_positive_data_set_file>

import sys
import os
from random import randint

if len(sys.argv) != 2:
	print "Usage: python random_positive_dataset.py <literature_positive_data_set_file>"
	print "Exiting..."
	sys.exit(1)

positive_data = open(sys.argv[1])

for line in positive_data:
	fields = line.split()
	print fields[randint(0, len(fields)-1)]

positive_data.close()
