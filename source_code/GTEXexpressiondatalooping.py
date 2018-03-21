#Requires GTEx expression data as input.

from math import log

with open("~/GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_metrics.tsv") as tissue_legend_tsv:
	tissue_legend_dict = {}
	tissue_legend_tsv.readline()
	current_line = tissue_legend_tsv.readline()
	while current_line:
		fields = current_line.split('\t')
		samplelabel = fields[0]
		tissue = fields[len(fields)-3]
		if tissue in tissue_legend_dict:
			tissue_legend_dict[tissue].append(samplelabel)
		else:
			tissue_legend_dict[tissue] = [samplelabel]        	
		current_line = tissue_legend_tsv.readline()

all_tissue_legend_keys = tissue_legend_dict.keys()

#Getting the entry locations associated with each key

with open("~/GTEx_Data_2014-01-17_RNA-seq_Flux1.6_transcript_reads.txt") as gtex_expression_data:
	current_line = gtex_expression_data.readline()
	column_labels = current_line.split()
	#print len(column_labels)
	tissue_locations_dict = {}
for key in all_tissue_legend_keys:
	all_patient_samples = tissue_legend_dict[key]
	patient_sample_locations_tissue = []
	for patient_sample in all_patient_samples:
		patient_sample_location = column_labels.index(patient_sample)
		patient_sample_locations_tissue.append(patient_sample_location)
	tissue_locations_dict[key] = patient_sample_locations_tissue

#Averaging the data associated with each key:

with open("~/GTEx_Data_2014-01-17_RNA-seq_Flux1.6_transcript_reads.txt") as gtex_expression_data:
	#this dict structure will contain the data
	print '\t'.join(["TranscriptID", '\t'.join(map(str,all_tissue_legend_keys)), "Tissue Entropy"])
	#skip the first line
	gtex_expression_data.readline()
	#and the second line
	current_line = gtex_expression_data.readline()
	while current_line:
		tissue_means_list = []
		tissue_contributions_entropy = []
		if current_line.startswith('TargetID'):
			current_line = gtex_expression_data.readline()
		if not current_line:
			break
		for key in tissue_locations_dict:
			fields = current_line.split()
			transcript_name = fields[0]
			expression_float = [float(fields[x]) for x in tissue_locations_dict[key]]
			tissue_mean = (sum(expression_float)/len(expression_float))
			tissue_means_list.append(tissue_mean)
		try:
			relative_expression = [float(x/sum(tissue_means_list)) for x in tissue_means_list]
			for x in relative_expression:
				if x == 0:
					tissue_contribution_entropy = 0
				else:
					tissue_contribution_entropy = -x*log(x, 2)
				tissue_contributions_entropy.append(tissue_contribution_entropy)
			tissue_entropy = sum(tissue_contributions_entropy)
		except ZeroDivisionError:
			relative_expression = float('Inf')
			tissue_entropy = 0
		print '\t'.join([transcript_name, '\t'.join(map(str,tissue_means_list)), str(tissue_entropy)])
		current_line = gtex_expression_data.readline()
