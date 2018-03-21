### Requires the discretization library for R.

print(paste0("Began running program: ", Sys.time()))

args<-commandArgs(TRUE)
 
setwd("~/")

date = Sys.time()

CODON = args[1]

# Read the raw data table
uORF_table_file = 			paste("uorfs_stats_table_", CODON, "_reannealed.txt", sep='')
uORF_table_uorfs_file = 	paste("uorfs_stats_table_", CODON, "_uORFsonly.txt", sep='')
uORF_table_R_file = 		paste("uorfs_stats_table_", CODON, "_R.txt", sep='')
uORF_table_unique_file = 	paste("uorfs_stats_table_", CODON, "_unique.txt", sep='')
#Triples:
leeandfritschandgao.pos.set = 	paste("leeandfritschandgao_pos_set_", CODON, ".txt", sep='')
#Doubles:
leeandfritsch.pos.set = 	paste("leeandfritsch_pos_set_", CODON, ".txt", sep='')
leeandgao.pos.set = 		paste("leeandgao_pos_set_", CODON, ".txt", sep='')
fritschandgao.pos.set = 	paste("fritschandgao_pos_set_", CODON, ".txt", sep='')
#Positives:
doubsortrips.pos.set = 		paste("doubsortrips_pos_set_", CODON, ".txt", sep='')
#Singles:
fritschonly.pos.set = 		paste("fritsch_pos_set_", CODON, ".txt", sep='')
leeonly.pos.set = 			paste("lee_pos_set_", CODON, ".txt", sep='')
gaoonly.pos.set = 			paste("gao_pos_set_", CODON, ".txt", sep='')
#Neutrals:
leeorfritschorgao.pos.set = paste("leeorfritschorgao_pos_set_", CODON, ".txt", sep='')

ribprof.pos.set = 			paste("ribprof_pos_set_", CODON, ".txt", sep='')
if (CODON == 'ATG') literature.pos.set = 		paste("~/uORF_stats/Positive_", CODON, "_uORF_Dataset_IDs.txt", sep='')
rand.literature.pos.set =	paste("random_lit_pos_set_", CODON, ".txt", sep='')

script.path = 				"~/"
stats.table.path = 			"~/uORF_stats/"
uORFs.seq.path = 			"~/uORFs/"

# Here mix and match alternative start codon datasets as necessary, maintaining the above names
# for correctness with my downstream code

# Processing out of n-terminal extensions
system(paste("python", paste(script.path, "uorf_nterm_filter.py -u", sep=''), 
			paste(stats.table.path, uORF_table_file, sep=''), 
			paste(uORFs.seq.path, "uORFs_", CODON, ".fa >", sep=''), 
			uORF_table_uorfs_file))

# Processing table for analysis in R
system(paste("python", paste(script.path, "uorfs_statstable_to_Rtable.py", sep=''),
			uORF_table_uorfs_file, ">", 
			uORF_table_R_file))

#Triples, Choose out those uORFs expressed in both Lee AND Fritsch AND Gao:
system(paste('grep "yes[[:space:]]yes[[:space:]]yes[[:space:]]no[[:space:]]no[[:space:]]no"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", leeandfritschandgao.pos.set))

#Doubles:

system(paste('grep "yes[[:space:]]yes[[:space:]]no[[:space:]]no[[:space:]]no[[:space:]]N_A"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", leeandfritsch.pos.set))
system(paste('grep "yes[[:space:]]no[[:space:]]yes[[:space:]]no[[:space:]]N_A[[:space:]]no"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", fritschandgao.pos.set))
system(paste('grep "no[[:space:]]yes[[:space:]]yes[[:space:]]N_A[[:space:]]no[[:space:]]no"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", leeandgao.pos.set))

system(paste('cat', leeandfritschandgao.pos.set, leeandfritsch.pos.set, fritschandgao.pos.set, leeandgao.pos.set, '>', doubsortrips.pos.set, sep=' '))

#Singles

system(paste('grep "yes[[:space:]]no[[:space:]]no[[:space:]]no[[:space:]]N_A[[:space:]]N_A"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", fritschonly.pos.set))
system(paste('grep "no[[:space:]]yes[[:space:]]no[[:space:]]N_A[[:space:]]no[[:space:]]N_A"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", leeonly.pos.set))
system(paste('grep "no[[:space:]]no[[:space:]]yes[[:space:]]N_A[[:space:]]N_A[[:space:]]no"', paste(stats.table.path, uORF_table_file, sep=''),
			 "| cut -f 1 >", gaoonly.pos.set))

# Assemble all of the positive examples, in the same file

system(paste('cat', leeonly.pos.set, fritschonly.pos.set, gaoonly.pos.set, '>', leeorfritschorgao.pos.set, sep=' '))

# Choose out a unique positive data set from literature uORF experiments
if (CODON == 'ATG') system(paste('python', paste(script.path, "random_positive_dataset.py", sep=''),
			literature.pos.set, ">", 
			rand.literature.pos.set))
		
# 18 February 2014 modification: Let downstream uniquing remove duplicates from this set.
if (CODON == 'ATG') system(paste('python', paste(script.path, "all_positive_dataset.py", sep=''),
			literature.pos.set, ">", 
			rand.literature.pos.set))

V = 1
TABLE_POSITIVE=paste("uorfs_stats_table_", CODON, "_R.positive", sep='')
TABLE_UNLABELED=paste("uorfs_stats_table_", CODON, "_R.unlabeled", sep='')
TABLE_NEUTRAL=paste("uorfs_stats_table_", CODON, "_R.neutral", sep='')

# Make positive and unlabeled datasets from the literature uORF experiments and ribosome profiling data
if (CODON == 'ATG') { system(paste('python', paste(script.path, "make_positive_and_negative_set_experiment.py", sep=''),
			"-d", uORF_table_R_file,
			"-r", leeorfritschorgao.pos.set,
			"-p", doubsortrips.pos.set,
			"-p", rand.literature.pos.set))
} else { system(paste('python', paste(script.path, "make_positive_and_negative_set_experiment.py", sep=''),
			"-d", uORF_table_R_file,
			"-r", leeorfritschorgao.pos.set,
			"-p", doubsortrips.pos.set))	
}

# Processing to get unique uORF tables without uORF start sites overlapping with CDSs
# This random choice is not exactly reproducible unless I set a seed in my python code, which I won't do for now.

INTERSECT="~/bin/intersectBed"
uORF_BED=paste("uORFs_", CODON, "_start.bed", sep='')
CDS_BED=paste("CDSs", CODON, "iteration.bed", sep='')
CDS_GTF="~/CDS/CDSs.gtf"

uORF_GTF=paste("~/uORFs/uORFs_", CODON, ".gtf", sep='')
INTERSECT_FILE=paste("uORFs_", CODON, "_CDS.intersect", sep='')
NONOVERLAP_CDS_uORF_IDs=paste("uORFs_", CODON, "_nonoverlap_ids.txt", sep='')
UNIQUE_IDs=paste("uORFs_", CODON, "_uniqids_", V, ".txt", sep='')
TMP_POS=paste("filtered_table_", CODON, "_positive.txt", sep='')
TMP_NEU=paste("filtered_table_", CODON, "_neutral.txt", sep='')
TMP_UNL=paste("filtered_table_", CODON, "_unlabeled.txt", sep='')
NUM_LIST=paste("~/Sandbox/numberedlist_", CODON, ".txt", sep='')
TMP=paste("tempdir_", CODON, ".txt", sep='')
TMP_NUM_POS=paste("filtered_numbered_table_", CODON, "_positive.txt", sep='')
TMP_NUM_NEU=paste("filtered_numbered_table_", CODON, "_neutral.txt", sep='')
TMP_NUM_UNL=paste("filtered_numbered_table_", CODON, "_unlabeled.txt", sep='')
FINAL_TABLE_POS=paste("filtered_table_positive_unique_", CODON, "_final.txt", sep='')
FINAL_TABLE_UNL=paste("filtered_table_unlabeled_unique_", CODON, "_final.txt", sep='')
FINAL_TABLE_NEU=paste("filtered_table_neutral_unique_", CODON, "_final.txt", sep='')

system(paste("python ~/uorf_gtf_to_BED_start.py", CODON, uORF_GTF, "."))
system(paste("cat ", CDS_GTF, " | awk -F $'\t' 'BEGIN {OFS = FS} ; {print $1,$4,$5,\"name\",$6,$7}' > ", CDS_BED))
system(paste(INTERSECT, " -a ", uORF_BED, " -b ", CDS_BED, " -wao > ", INTERSECT_FILE))

print(paste0("Completed intersect of CDS and uORF start: ", Sys.time()))

# #This next program, because of the sys.stderr.write command, prints a bunch of uORF IDs at the command line. Not helpful visually.

system(paste("python ~/seg_uORFids_CDSoverlap.py ", INTERSECT_FILE, " > ", NONOVERLAP_CDS_uORF_IDs))
system(paste("python ~/rand_uniq_uORF_IDs_pos_neu.py ", uORF_GTF, " > ", UNIQUE_IDs))
system(paste("python ~/filter_uORFs.py ", TABLE_POSITIVE, NONOVERLAP_CDS_uORF_IDs, " > ", TMP_POS))
system(paste("python ~/filter_uORFs.py", TMP_POS, UNIQUE_IDs, " > ", FINAL_TABLE_POS))

system(paste("python ~/seg_uORFids_CDSoverlap.py ", INTERSECT_FILE, " > ", NONOVERLAP_CDS_uORF_IDs))
system(paste("python ~/rand_uniq_uORF_IDs.py ", uORF_GTF, " > ", UNIQUE_IDs))
system(paste("python ~/filter_uORFs.py ", TABLE_UNLABELED, " ", NONOVERLAP_CDS_uORF_IDs, " >", TMP_UNL))
system(paste("python ~/filter_uORFs.py", TMP_UNL, UNIQUE_IDs, " > ", FINAL_TABLE_UNL))

system(paste("python ~/seg_uORFids_CDSoverlap.py ", INTERSECT_FILE, " > ", NONOVERLAP_CDS_uORF_IDs))
system(paste("python ~/rand_uniq_uORF_IDs.py ", uORF_GTF, " > ", UNIQUE_IDs))
system(paste("python ~/filter_uORFs.py ", TABLE_NEUTRAL, " ", NONOVERLAP_CDS_uORF_IDs, " >", TMP_NEU))
system(paste("python ~/filter_uORFs.py", TMP_NEU, UNIQUE_IDs, " > ", FINAL_TABLE_NEU))

#*****************below should not be commented

system(paste("python ~/seg_uORFids_CDSoverlap.py ", INTERSECT_FILE, " > ", NONOVERLAP_CDS_uORF_IDs))
system(paste("python ~/rand_uniq_uORF_IDs_pos_neu.py ", uORF_GTF, " > ", NUM_LIST))

TMP_POS=paste("filtered_table_positive_", CODON, ".txt", sep='')
TMP_NEU=paste("filtered_table_neutral_", CODON, ".txt", sep='')
TMP_UNL=paste("filtered_table_unlabeled_", CODON, ".txt", sep='')

system(paste("python ~/filter_uORFs.py ", TABLE_POSITIVE, " ", NONOVERLAP_CDS_uORF_IDs, " > ", TMP_POS))
system(paste("python ~/filter_uORFs.py ", TABLE_NEUTRAL, " ", NONOVERLAP_CDS_uORF_IDs, " > ", TMP_NEU))
system(paste("python ~/filter_uORFs.py ", TABLE_UNLABELED, " ", NONOVERLAP_CDS_uORF_IDs, " > ", TMP_UNL))

system(paste("python ~/filter_uORFs_2.py ", TMP_POS, " ", NUM_LIST, " > ", TMP_NUM_POS))
system(paste("python ~/filter_uORFs_2.py ", TMP_NEU, " ", NUM_LIST, " > ", TMP_NUM_NEU))
system(paste("python ~/filter_uORFs_2.py ", TMP_UNL, " ", NUM_LIST, " > ", TMP_NUM_UNL))

#The following is a uniquing step. The output of this step is three separate files.

system(paste("python ~/unique_pos_neu_unl_speed.py", CODON, sep=' '))

# The following pulls the data and moves it to final_table_unl.

system(paste("python ~/retrieveuniquedlines.py", CODON, sep=' '))

# Add the uORF table header back to the positive, unlabeled, and neutral data files.
system(paste("head -n +1", uORF_table_R_file, ">", TMP))
system(paste("cat", TMP, FINAL_TABLE_POS, "> tmp"))
system(paste("cat tmp >", FINAL_TABLE_POS))

system(paste("head -n +1", uORF_table_R_file, ">", TMP))
system(paste("cat", TMP, FINAL_TABLE_UNL, "> tmp"))
system(paste("cat tmp >", FINAL_TABLE_UNL))

system(paste("head -n +1", uORF_table_R_file, ">", TMP))
system(paste("cat", TMP, FINAL_TABLE_NEU, "> tmp"))
system(paste("cat tmp >", FINAL_TABLE_NEU))

# system(paste("rm", TMP_POS))
# system(paste("rm", TMP_UNL))
# system(paste("rm", TMP_NEU))
system(paste("rm", TMP))
system(paste("rm", "tmp"))

print(paste0("Removed temporary files: ", Sys.time()))

# Load the fully processed data tables into R
print(paste0("Loading the positive data table: ", Sys.time()))
positive_data =  read.table(FINAL_TABLE_POS, header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
print(paste0("Loading the unlabeled data table: ", Sys.time()))
unlabeled_data = read.table(FINAL_TABLE_UNL, header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
print(paste0("Loading the neutral data table: ", Sys.time()))
neutral_data = read.table(FINAL_TABLE_NEU, header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")

# Save the positive data set used for future reference, in repeating this procedure
write.table(positive_data, file=paste("processed_positiveset_", CODON, "_uORF_table_", date, ".txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE)

# Combine the positive dataset and unlabeled datasets and discretize them

print("Reached Discretization step")
library(discretization, lib.loc="~/")
setwd("~/Rfiles")

# The mdlp discretization algorithm uses the last column of the matrix as the class variables
pos.d = cbind(positive_data, rep(1, dim(positive_data)[1]))
unl.d = cbind(unlabeled_data, rep(0, dim(unlabeled_data)[1]))
neu.d = cbind(neutral_data, rep(2, dim(neutral_data)[1]))
colnames(pos.d)[dim(pos.d)[2]] <- "class.label"
colnames(unl.d)[dim(unl.d)[2]] <- "class.label"
colnames(neu.d)[dim(unl.d)[2]] <- "class.label"
all.data = rbind(pos.d, unl.d, neu.d)

print(paste0("Running mdlp algorithm: ", Sys.time()))

all.disc = mdlp(all.data)$Disc.data

print(paste0("Discretized the data: ", Sys.time()))

uORF.c = all.data["uORF.context"] + 1			# Add one to make the values work with Naive Bayes code
c.lab = all.data["class.label"]

all.disc = subset(all.disc, select = -uORF.context) 
all.disc = cbind(all.disc, uORF.c)
colnames(all.disc)[dim(all.disc)[2]] = "uORF.context"
all.disc = subset(all.disc, select = -class.label)
all.disc = cbind(all.disc, c.lab)
colnames(all.disc)[dim(all.disc)[2]] = "class.label"

print(paste0("Restored uORF context variable: ", Sys.time()))

write.table(all.disc, file=paste("discretized_", CODON, "_uORF_table_", date, sep=""), quote=FALSE, sep="\t", col.names=TRUE)

print(paste0("Wrote the discretized data in a file: ", Sys.time()))
