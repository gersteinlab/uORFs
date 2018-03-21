import os
import sys

specific_codon = sys.argv[1]

output_file = open("~/uORF_stats/uorfs_stats_table_" + specific_codon + "_reannealed.txt", 'w')

filenames = ["~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition1.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition2.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition3.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition4.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition5.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition6.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition7.txt", \
"~/uORF_stats/uorfs_stats_table_" + specific_codon + "_partition8.txt"]

header_file = open(filenames[0], 'r')
header = header_file.readline()
output_file.write(header) 

for fname in filenames:
	with open(fname) as infile:
		for line in infile:
				if not line.startswith('#'):
					output_file.write(line)
