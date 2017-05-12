import os

all_start_codons = open("~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/allcodons_inuorf.txt", 'r')

codon_list = []
for item in all_start_codons:
	codon = item.split()[0].strip()
	codon_list.append(codon)
print codon_list

uorfs_script = "python ~/python_scripts/uorfs.py"
annotation_file = "/net/gerstein/GENCODE/gencode19/gencode.v19.annotation.gtf"
chromsomes_dir = "/net/gerstein/genomes/human/hg19/igenome/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/"
output_directory = "~/uORFs/"

uorf_stats_bashscript = "~/bash_scripts/uorf_stats_table.bash"
lee_IDs = "~/ribosome_profiling/globalmapping_uORFs_paper_Leeetal2012/lee_ATG_uORF_IDs"
fritsch_IDs = "~/ribosome_profiling/novelhuman_uORFs_paper_Fritschetal2012/UCSCBrowser_identifiedTISs/fritsch_ATG_uORF_IDs"

# # for codon in codon_list:
# # 	os.system(uorfs_script + ' ' + annotation_file + ' ' + chromsomes_dir + ' ' + output_directory + ' ' + '-' + codon)

# gtf_to_bed_script = "python ~/python_scripts/uorf_gtf_to_bed.py"

# for codon in codon_list:
# 	os.system(gtf_to_bed_script + ' ' + codon)

# SNP_intersections_script = "python ~/python_scripts/SNPdataextract.py"

# for codon in codon_list:
# 	os.system(SNP_intersections_script + ' ' + codon)

# for codon in codon_list:
os.system('sh ' + uorf_stats_bashscript + ' ' + 'AAG' + ' ' + lee_IDs + ' ' + fritsch_IDs)

# for codon in codon_list:
# 	os.system("python ~/python_scripts/number_uORFs_per_start_codon.py" + ' ' + codon)