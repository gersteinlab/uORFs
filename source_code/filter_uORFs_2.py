import os
import sys

# system(paste("python ~/filter_uORFs.py", TMP_POS, UNIQUE_IDs, " > ", FINAL_TABLE_POS))

if __name__=='__main__':

	uORF_IDs = {}
	#open the uORF data file
	uORF_data_file = open(sys.argv[1])
	#perform the add_filter_IDs function on UNIQUE_IDs. This basically adds all the IDs to a dict.
	file = open(sys.argv[2])
	for line in file:
		uORF_IDs[line.split()[1]] = line.split()[0]
	for line in uORF_data_file:
		if line.split()[0] in uORF_IDs:
				print '\t'.join([line.split()[0], uORF_IDs[line.split()[0]]])
		else:
			if line.split()[0] not in uORF_IDs:
				print '\t'.join([line.split()[0], 'N'])
	uORF_data_file.close()
	 
		
