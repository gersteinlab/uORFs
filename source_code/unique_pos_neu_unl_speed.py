import os
import sys
from random import randint

if __name__=='__main__':

	CODON = sys.argv[1]

	positive_read_file = "~/Sandbox/filtered_numbered_table_" + CODON + "_positive.txt"
	positive_write_file = "~/Sandbox/filtered_numbered_table_positive_" + CODON + "_unique.txt"
	neutral_read_file = "~/Sandbox/filtered_numbered_table_" + CODON + "_neutral.txt"
	neutral_write_file = "~/Sandbox/filtered_numbered_table_neutral_" + CODON + "_unique.txt"
	unlabeled_read_file = "~/Sandbox/filtered_numbered_table_" + CODON + "_unlabeled.txt"
	unlabeled_write_file = "~/Sandbox/filtered_numbered_table_unlabeled_" + CODON + "_unique.txt"

	fo = open(positive_write_file, "w")
	uORF_data_file = open(positive_read_file)
	already_covered = set()
	for line in uORF_data_file:
		number = line.split('\t')[1]
		number = number.strip()
		if number not in already_covered:
			already_covered.add(number)
			randcollection = list()
			uORF_data_filetwo = open(positive_read_file)
			for entry in uORF_data_filetwo:
				numbercompare = entry.split('\t')[1]
				numbercompare = numbercompare.strip()
				ID = entry.split('\t')[0]
				ID = ID.strip()
				if numbercompare == number:
					randcollection.append(ID)
			if len(randcollection) >= 1:
				fo.write(''.join([randcollection[randint(0, len(randcollection)-1)], '\n']))
			else:
				fo.write(''.join([randcollection[0], '\n']))
	
	fo = open(neutral_write_file, "w")
	uORF_data_file = open(neutral_read_file)
	for line in uORF_data_file:
		number = line.split('\t')[1]
		number = number.strip()
		if number not in already_covered:
			already_covered.add(number)
			randcollection = list()
			uORF_data_filetwo = open(neutral_read_file)
			for entry in uORF_data_filetwo:
				numbercompare = entry.split('\t')[1]
				numbercompare = numbercompare.strip()
				ID = entry.split('\t')[0]
				ID = ID.strip()
				if numbercompare == number:
					randcollection.append(ID)
			if len(randcollection) >= 1:
				fo.write(''.join([randcollection[randint(0, len(randcollection)-1)], '\n']))
			else:
				fo.write(''.join([randcollection[0], '\n']))

	fo = open(unlabeled_write_file, "w")
	uORF_data_file = open(unlabeled_read_file)
	for line in uORF_data_file:
		number = line.split('\t')[1]
		number = number.strip()
		if number not in already_covered:
			already_covered.add(number)
			randcollection = list()
			uORF_data_filetwo = open(unlabeled_read_file)
			for entry in uORF_data_filetwo:
				numbercompare = entry.split('\t')[1]
				numbercompare = numbercompare.strip()
				ID = entry.split('\t')[0]
				ID = ID.strip()
				if numbercompare == number:
					randcollection.append(ID)
			if len(randcollection) >= 1:
				fo.write(''.join([randcollection[randint(0, len(randcollection)-1)], '\n']))
			else:
				fo.write(''.join([randcollection[0], '\n']))
