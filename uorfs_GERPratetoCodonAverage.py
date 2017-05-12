#! /usr/bin/env python

# Author: Russell Ault

# uorfs_GERPratetoCodonAverage.py <GERPrate_start_stop file>

# The program output goes to the standard output
# Output:
# uORF_ID	start	end
# ENST099	1.333	2.3333

# Gather each ID's start and stop lines into two separate lists

import sys

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print "Incorrect Number of Arguments"
		print "Usage: uorfs_GERPratetoCodonAverage.py <GERPrate_start_stop file>"
		print "Exiting..."
		sys.exit(1)
	
	# Header Line
	print "#uORF_ID\tstartCodon-GERP-score-average\tstopCodon-GERP-score-average"
	
	GERPrateInputFileName = sys.argv[1]
	with open(GERPrateInputFileName) as GERPrateInputFile:
		prevID = None
		start_scores = []
		stop_scores = []
		line=GERPrateInputFile.readline()
		while line:
			fields = line.split()
			ID = fields[4]
			if ID != prevID and prevID:
				print '\t'.join([prevID, str(sum(start_scores)/float(len(start_scores))),
									str(sum(stop_scores)/float(len(stop_scores)))])
				start_scores = []
				stop_scores = []
			

			if 'start' in line:
				start_scores.append(float(fields[6]))
			elif 'stop' in line:
				stop_scores.append(float(fields[6]))
			else:
				print "Incorrect input line format error, no start or stop string in line"
				print "Exiting..."
				sys.exit(1)
			prevID = ID
			line=GERPrateInputFile.readline()
			
			if not line:					# Processes last uORF
				print '\t'.join([prevID, str(sum(start_scores)/float(len(start_scores))),
									str(sum(stop_scores)/float(len(stop_scores)))])
				break
	
									
			