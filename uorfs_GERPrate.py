#! /usr/bin/env python

# Author: Russell Ault

# uorfs_GERPrate.py <somekindofgerprate directiory> <some kind of gerp cache directory> -v <VATty file format for my uORFstarts and stops> -i <interval file> <-- This last one will not be used yet

# uorfs_GERPrate.py takes a file in .gerp format (similar to the first few lines of VAT) and calculates the GERP
# score average for the bases represented by each line. Further processing is required to average start site and stop site
# bases for individual genomic uORFs.

import os
import sys
from common import *

global gerpOutputFile
global Codon

# WRITE THIS LATER
def getIntervalData(chrs, GERPelementpath, codingExonIntervals):
	pass

def getGERPData(gerpInputFile, chrs, GERPelementpath, GERPratepath, GERPratecachepath, codingExonIntervals=None):
	## Coordinates are 1-based.
	## All GERP intervals include endpoints
	
	gerpInputFile.seek(0)
	line = gerpInputFile.readline()
	while line.startswith("#") or line=="\n":
		line=gerpInputFile.readline()					# skip any intro lines or blank lines

	for i in chrs:
		gerpCacheFile = buildGerpRates(GERPratepath, GERPratecachepath, i)
		while line.split('\t')[0].split('chr')[-1]==i:
			data = line.split()
			#chr_num = data[0].split('chr')[-1]
			chromosome = data[0]
			start = int(data[1])
			length = int(data[2])
			end = start + length-1  ##inclusive endpoint
			strand = data[3]
			uORF_ID  = data[4]
			start_or_stop = data[5]
			
			GERPrateValue = getGerpScore(gerpCacheFile, start, length)
			gerpOutputFile.write('\t'.join([chromosome, str(start), str(length), strand, uORF_ID, 
							  start_or_stop, str(GERPrateValue)]) + '\n')
			line=gerpInputFile.readline()

if __name__== '__main__':

#uorfs_GERPrate.py CODON <somekindofgerprate directiory> <some kind of gerp cache directory> -v <VATty file format for my uORFstarts and stops> -i <interval file>
	if len(sys.argv) != 6:
		print "Incorrect Number of Arguments!"
		print "Usage: uorfs_GERPrate.py <somekindofgerprate directiory> <some kind of gerp cache directory> [-v <VATty file format for my uORFstarts and stops>, -i <interval file>]"
		print "Exiting..."
		sys.exit(1)
		
	chromosome_lst = [str(i) for i in range(1,23)] + ['X', 'Y']
	
	GERPelementpath = None
	Codon = sys.argv[1]
	GERPratePath = sys.argv[2]
	GERPcachePath = sys.argv[3]
	if sys.argv[4] == "-v":
		gerpInputFileName = sys.argv[5]
		gerpOutputFile = open("~/GERP_results/GERPrate_start_stop_" + Codon, "w")
		with open(gerpInputFileName) as gerpInputFile:
			getGERPData(gerpInputFile, chromosome_lst, GERPelementpath, GERPratePath, GERPcachePath)
		gerpOutputFile.close()
			
	elif sys.argv[4] == "-i":
		intervalInputFileName = sys.argv[5]
		with open(intervalInputFileName) as intervalInputFile:
			getIntervalData() # Fill this in later!