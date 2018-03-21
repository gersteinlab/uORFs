### Requires 1KG variant vcf as input, BEDTools intersect as dependency

import gzip
import sys
import os

if __name__ == '__main__':

	codon = sys.argv[1]

	tracker = 1
	SNP_data = gzip.open('~/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz', 'r')
	SNP_dict = {}
	#skip the first line
	SNP_data.readline()
	#store the second line
	current_line = SNP_data.readline()
	f = open("~/python_scripts/bedresultSNP.bed",'w')
	while current_line:
		if current_line.startswith('#'):
			current_line = SNP_data.readline()
		else:
			break
	while current_line:
		fields = current_line.split()
		SNPstatsfields = fields[7].split(";AF=")[1]
		AFscore = SNPstatsfields.split(';')[0]
		chromo = str(fields[0])
		chromobed = ''.join(["chr", chromo])
		posit = fields[1]
		positminusone = int(posit)-1
		namefield = fields[1]
		# patient_counts_floats = [float(x) for x in patient_counts]
		# patient_counts_floats_mean = (sum(patient_counts_floats)/len(patient_counts_floats))
		# gtex_expression_dict[transcript_name] = patient_counts_floats_mean
		# current_line = gtex_expression_data.readline()
		chromo = fields[0]
		posit = fields[1]
		print >> f, '\t'.join([chromobed, str(positminusone), str(posit), str(AFscore)])
		current_line = SNP_data.readline()
	#		tracker += 1
	#	else:
	f.close()

	next line, sorts the file, for improved speed, if necessary.
	os.system("sort -k1,1 -k2,2n ~/uORFs/g19/uORF_ATG.bed > ~/python_scripts/uORF_ATG_sorted.bed")
# 
	os.system('STARTTIME=$(date) && echo $STARTTIME && STARTTIMESECS=$(date +%s)')
	os.system('~/bin/intersectBed -wo -a ~/uORFs/uORF_' + codon + '.bed -b ~/python_scripts/bedresultSNP.bed > ~/SNPs_intersections/' + codon + '_SNPs_uORFs.intersect')

	tracker = 1
	with open("~/SNPs_intersections/" + codon + "_SNPs_uORFs.intersect") as intersect_data:
		#this dict structure will contain the data
		intersect_dict = {}
		# read the first line
		current_line = intersect_data.readline()
		while current_line:
			if current_line == '\n':
				break
			if not current_line:
				break
			fields = current_line.split()
			transcript_ID = fields[3]
			mutation_frequency = fields [-2]
			if transcript_ID in intersect_dict:
				intersect_dict[transcript_ID].append(mutation_frequency)
				current_line = intersect_data.readline()
			else:
				intersect_dict[transcript_ID] = [mutation_frequency]
				current_line = intersect_data.readline()
			tracker += 1

	hetero_snp_write_file = open("~/SNPs_intersections/" + codon + "_heteroSNPfile.txt", 'w')
	with open("~/uORFs/uORFs_" + codon + ".fa") as uORFs_seq:
		uORFs_seq.seek(0)
		uORFs_seq.readline()		# skip intro line
		data = uORFs_seq.readline()
		print >> hetero_snp_write_file, "\t".join(["uORFID", "NPlyMrphsm", "Htrzygty"])
		while data:
			fields = data.split('|')
			uORF_ID = fields[2]
			uORFs_seq.readline() # Skip the sequence line.
			mutation_frequencies = intersect_dict.get(uORF_ID)
			if mutation_frequencies is not None:
				heteroscore = 0
				heterosum = 0
				wildtype = 0
				howmany = 0
				for x in mutation_frequencies:
					heteroscore = heteroscore + (float(x))**2 + (1-float(x))**2
					howmany += 1
				heteroscore = 1-((heteroscore)/howmany)
				print >> hetero_snp_write_file, "\t".join([uORF_ID, str(howmany), str(heteroscore)])
			else:
				print >> hetero_snp_write_file, "\t".join([uORF_ID, "0", "0"])
			data = uORFs_seq.readline()
			tracker += 1
			#Skip header lines.
			if data and data.startswith('#'):
				data = uORFs_seq.readline()
			tracker += 1
