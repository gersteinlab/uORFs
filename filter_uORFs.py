# python filter_uORFs.py [-r] [-p ID-POSITION] <uORF-data-file> [uORF-ID-file]*

# take into account a comment line

import os
import sys
from make_negativeset import ERROR


def add_filter_IDs(filename, IDs=set(), ID_pos=0):
	file = open(filename)
	for line in file:
		IDs.add(line.split()[ID_pos])
	return IDs


if __name__=='__main__':
	if len(sys.argv) < 3:
		ERROR("Usage: python filter_uORFs.py [-r] [-p ID-POSITION] <uORF-data-file> [uORF-ID-file]*")
	
	i = 1
	ID_pos = 0
	uORF_IDs_to_filter = set()
	keep_IDs = True
	if sys.argv[i] == '-r':
		keep_IDs = False
		i += 1
	if sys.argv[i] == '-p':
		ID_pos = int(sys.argv[i+1])
		i += 2
	uORF_data_file = open(sys.argv[i])
	i += 1
	while i < len(sys.argv):
		uORF_IDs_to_filter = add_filter_IDs(sys.argv[i], uORF_IDs_to_filter)
		i += 1
	for line in uORF_data_file:
		if not line.startswith('#'):
			if keep_IDs:
				if line.split()[ID_pos] in uORF_IDs_to_filter:
					print line,
			else:
				if line.split()[ID_pos] not in uORF_IDs_to_filter:
					print line,
	uORF_data_file.close()
	 
		