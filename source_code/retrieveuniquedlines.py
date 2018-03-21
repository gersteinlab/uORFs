import os
import sys

if __name__=='__main__':

	CODON = sys.argv[1]

#Retrieval and writing of positive examples:

	positive_filtered_file = "filtered_table_" + CODON + "_positive.txt"
	positive_unique_file =  "~/Sandbox/filtered_numbered_table_positive_" + CODON + "_unique.txt"
	positive_unique_final_file =  "~/Sandbox/filtered_table_positive_unique_" + CODON + "_final.txt"
	neutral_filtered_file = "filtered_table_" + CODON + "_neutral.txt"
	neutral_unique_file =  "~/Sandbox/filtered_numbered_table_neutral_" + CODON + "_unique.txt"
	neutral_unique_final_file =  "~/Sandbox/filtered_table_neutral_unique_" + CODON + "_final.txt"
	unlabeled_filtered_file = "filtered_table_" + CODON + "_unlabeled.txt"
	unlabeled_unique_file =  "~/Sandbox/filtered_numbered_table_unlabeled_" + CODON + "_unique.txt"
	unlabeled_unique_final_file =  "~/Sandbox/filtered_table_unlabeled_unique_" + CODON + "_final.txt"


	all_IDs = open(positive_filtered_file, "r")
	with open('~/Sandbox/uorfs_stats_table_ATG_uORFsonly.txt', 'r') as f:
    		headerline = f.readline()
	
	output_file = open(positive_unique_final_file, "w")
	output_file.write(headerline)
	for line in all_IDs:
		ID = line.split('\t')[0]
		unique_IDs = open(positive_unique_file, "r")
		for item in unique_IDs:
			item = item.strip()
			if item == ID:
				output_file.write(line)

#Retrieval and writing of unlabeled examples:

	all_IDs = open(unlabeled_filtered_file, "r")
	output_file = open(unlabeled_unique_final_file, "w")
	output_file.write(headerline)
	for line in all_IDs:
		ID = line.split('\t')[0]
		unique_IDs = open(unlabeled_unique_file, "r")
		for item in unique_IDs:
			item = item.strip()
			if item == ID:
				output_file.write(line)

# Retrieval and writing of neutral examples:

	all_IDs = open(neutral_filtered_file, "r")
	output_file = open(neutral_unique_final_file, "w")
	output_file.write(headerline)
	for line in all_IDs:
		ID = line.split('\t')[0]
		unique_IDs = open(neutral_unique_file, "r")
		for item in unique_IDs:
			item = item.strip()
			if item == ID:
				output_file.write(line)
