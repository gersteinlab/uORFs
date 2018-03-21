# python rand_uniq_uORF_IDs.py <uORFs.gtf file>

# rand_uniq_uORF_IDs.py takes a uORF GTF file and outputs a random selection of uORF IDs
# that are all unique (not the same chromosomal coordinates as any other uORF)

# This file also depends on the existence of the uorf_gtf_to_PhyloCSF.py script
# Output is to standard output.

import os
import sys

from make_negativeset import ERROR
from random import randint


def get_intervals(coord):
	intervals = [section.split(':')[1] for section in coord.split('+')]
	return tuple([(interval.split('-')[0], interval.split('-')[1]) for interval in intervals])
					

def all_uORFs_from_Phylo(PhyloCSF_file, coord_pos=0, strand_pos=1, ID_pos=2):
	
	dict = {}
	
	for line in PhyloCSF_file:
		fields = line.split()
		chromosome = fields[coord_pos].split(':')[0]
		strand = fields[strand_pos]
		uORF_ID = fields[ID_pos]
		intervals = get_intervals(fields[coord_pos])
		key = (chromosome, strand, intervals)
		
		if dict.get(key):
			IDs = list(dict[key])
			IDs.append(uORF_ID)
			dict[key] = tuple(IDs)
		else:
			dict[key] = tuple([uORF_ID])
			
	PhyloCSF_file.close()
	return dict

""" Structure of the key:
		(chromosome,
		 strand,
		((block1_start,block2_end),
		 (block2_start,block2_end),
		 .
		 .
		 .))
		 
		Structure of value:
		(uORF_ID1,
		 uORF_ID2,
		 uORF_ID3)"""

def random_ID(IDs_lst):
	return IDs_lst[randint(0, len(IDs_lst)-1)]

	
if __name__=='__main__':
	if len(sys.argv) != 2:
		ERROR("Usage: python rand_uniq_uORF_IDs <uORFs.gtf file>")
	
	uORF_GTF_filename = sys.argv[1]
	
	# This dictionary contains the coordinates for each uORF in a tuple of tuples, along 
	# with the strand, as the key, and the value is a tuple of all of the uORF IDs associated with that unique uORF
	uORF_dict = all_uORFs_from_Phylo(os.popen(' '.join(["python", "~/uorf_gtf_to_PhyloCSF.py", uORF_GTF_filename])))

	i = 0
	for IDs_lst in uORF_dict.values():
		j = 0
		for j in range(0, len(IDs_lst)):
			print '\t'.join([str(i), IDs_lst[j]])
			j +=1
		i += 1
