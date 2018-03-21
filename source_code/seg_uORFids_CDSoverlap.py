# Author: Russell Ault

# seq_uORFids_CDSoverlap.y takes the result of applying intersectBed to uORF start sites
# and CDS exons of the GENCODE annotation and outputs non-overlapping uORF IDs to standard output
# and overlapping uORF IDs to standard error

# python seg_uORFids_CDSoverlap.py <uORFs_CDS.intersect file>

import os
import sys


def initial_dict(intersect_file, ID_pos=3):
	intersect_file.seek(0)
	dict = {}
	for line in intersect_file:
		dict[line.split()[ID_pos]] = 0
	return dict


if __name__=='__main__':
	if len(sys.argv) != 2:
		print "Usage: python seg_uORFids_CDSoverlap.py <uORFs_CDS.intersect file>"
		sys.exit(1)
	
	intersect_file = open(sys.argv[1])
	
	non_overlap_IDs = set()
	overlap_IDs = set()
	
	ID_pos = 3
	uORF_over_dict = initial_dict(intersect_file)
	intersect_file.seek(0)
	for line in intersect_file:
		fields = line.split()
		uORF_over_dict[fields[ID_pos]] += int(fields[-1])
	for ID, over_bases in uORF_over_dict.items():
		if over_bases == 0:
			print ID
		else:
			sys.stderr.write(ID + '\n')
	intersect_file.close()
	
	
"""	for line in intersect_file:
		fields = line.split()
		if int(fields[-1]) == 0:
			non_overlap_IDs.add(fields[ID_pos])
		else:
			overlap_IDs.add(fields[ID_pos])
	
	for ID in non_overlap_IDs:
		print ID
	for ID in overlap_IDs:
		sys.stderr.write(ID + '\n')
	
	intersect_file.close()"""
