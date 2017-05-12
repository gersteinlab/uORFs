import os
import sys

# specific_codon is only used on occasion, when running the uorfs_stats_bashscript

specific_codon = sys.argv[1]
partition_number = sys.argv[2]

all_start_codons = open("~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/allcodons_inuorf.txt", 'r')

codon_list = []
for item in all_start_codons:
	codon = item.split()[0].strip()
	codon_list.append(codon)
print codon_list

uorfs_script = "python ~/uorfs.py"
annotation_file = "~/gencode.v19.annotation.gtf"
chromsomes_dir = "~/" ### Replace with directory of chromosomal fasta files (sequences)
output_directory = "~/uORFs/"

uorf_stats_bashscript = "~/uorf_stats_table_partition.bash"
lee_IDs = "~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/lee_" + specific_codon + "_uORF_IDs"
fritsch_IDs = "~/ribosome_profiling/novelhuman_uORFs_paper_Fritschetal2012/UCSCBrowser_identifiedTISs/fritsch_" + specific_codon + "_uORF_IDs"
gao_IDs = "~/ribosome_profiling/uORFs_paper_Gaoetal2014/gao_" + specific_codon + "_uORF_IDs"
simple_queue_file = open("~/partitioned_uORF_stats", 'w')

#for codon in codon_list:
#	i = 1
#	while i <= 8:
#		print >> simple_queue_file, "cd ~/python_scripts; python all_start_codons_uORFs_partition.py " + codon + " " + str(i)
#		i += 1

# # for codon in codon_list:
# # 	os.system(uorfs_script + ' ' + annotation_file + ' ' + chromsomes_dir + ' ' + output_directory + ' ' + '-' + codon)

# gtf_to_bed_script = "python ~/python_scripts/uorf_gtf_to_bed.py"

# for codon in codon_list:
# 	os.system(gtf_to_bed_script + ' ' + codon)

# SNP_intersections_script = "python ~/python_scripts/SNPdataextract.py"

# for codon in codon_list:
# 	os.system(SNP_intersections_script + ' ' + codon)

# for codon in codon_list:
os.system('sh ' + uorf_stats_bashscript + ' ' + specific_codon + ' ' + partition_number + ' ' + lee_IDs + ' ' + fritsch_IDs + ' ' + gao_IDs)

# sh uorf_stats_bashscript specific_codon partition_number lee_IDs fritsch_IDs gao_IDs)
# 	os.system("python ~/python_scripts/number_uORFs_per_start_codon.py" + ' ' + codon)