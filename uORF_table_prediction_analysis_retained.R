date = Sys.Date()

setwd("~/")

#-------------------------------------
# Function Definitions

percent_pos = function(training, pred, class, freq=F) {
	# If both are 1 in the same place, then we have percent of positives
	n_class = sum(training == class)
	n_pos = sum(((training == class) + (pred == 1)) == 2)
	if (freq) return(n_pos)
	else return(n_pos / n_class)
}

get.perf.measures.yao = function(pred, labels) {

	
	p_pos_pos = percent_pos(labels, pred, T)
	p_unl_pos = percent_pos(labels, pred, F)
	
	# Precision calculation
	pos_pos = percent_pos(labels, pred, T, freq=T)
	unl_pos = percent_pos(labels, pred, F, freq=T)
	prec = pos_pos / (pos_pos + unl_pos)
	recall = p_pos_pos			# These are the same thing
	F.stat = 2 * (prec*recall) / (prec + recall)
	return(c(p_pos_pos, p_unl_pos, prec, recall, F.stat))

	# Precision = true-positives / (true-positives + false-positives)
	# recall = true-positives / (true-positives + false-negatives)

	# F statistic is the geometric mean of Precision and Recall:
	# F = 2 * (precision * recall) / (precision + recall)
	
}

n.Bayes = function(all.data, u.prb, p.prb, prior) {
	M = dim(all.data)[1]
	N = dim(all.data)[2]
	data  = all.data[,-N]
	labels = all.data[,N]
	pred = logical(M)
	
	for (i in 1:M) {
		pred[i] = pred.Bayes(as.numeric(data[i,]), u.prb, p.prb, prior)
	}
	print(paste("Finished one prior: ", prior))
	return(get.perf.measures.yao(pred, labels))
	
}

n.Bayes.score = function(data, u.prb, p.prb, prior) {
	scores = data[,1]			# To retain the row names
	M = dim(data)[1]
	for (i in 1:M) {
		scores[i] = pred.Bayes(as.numeric(data[i,]), u.prb, p.prb, prior, ratio=T)
	}
	
	return(scores)
}

# p is the probability of positive
pred.Bayes = function(x, u.prb, p.prb, p, ratio=FALSE) {
	p.unl = log10((1-p))
	n = dim(u.prb)[1]
	for (j in 1:n) {
		p.unl = p.unl + log10(u.prb[j, x[j]])
		# deb = p.prb[j, x[j]]
		# print(paste("Deb:", deb, "j:", j, "x[j]", x[j]))
		# if (x[j] == 1) {
			# p.unl = p.unl * u.prb[j]
		# }
		# else if (x[j] == 2) {
			# p.unl = p.unl * (1 - u.prb[j])
		# }
		# else print("ERROR, feature value not discretized as expected")
		
	}
	# print(p.prb)
	# print(u.prb)
	p.pos = log10(p)
	for (j in 1:n) {
		p.pos = p.pos + log10(p.prb[j, x[j]])
		# deb = p.prb[j, x[j]]
		# print(paste("Deb:", deb, "j:", j, "x[j]", x[j]))
		# if (x[j] == 1) {
			# p.pos = p.pos * p.prb[j]
		# }
		# else if (x[j] == 2) {
			# p.pos = p.pos * (1 - p.prb[j])
		# }
		# else print("ERROR, feature value not discretized as expected")
		
	}
	# print(p.pos)
	# print(p.unl)
	if (ratio == TRUE) {
		return(p.pos - p.unl)
	}
	if (p.pos > p.unl) return(TRUE)
	else return(FALSE)
	
}

screen = function(M) {
	m = dim(M)[1]
	n = dim(M)[2]
	bad.c = c()
	for (k in 1:n) {
		fr = mean(M[,k] == 1)
		if (fr == 0.0 || fr == 1.00) {bad.c = c(bad.c, k)}
	}
	return(M[, -bad.c])
	
}
sample.mat = function(M, p) {
	n = dim(M)[1]
	m = dim(M)[2]
	f = 100000
	p.data = as.integer(p * f)
	log = c(rep(T, p.data), rep(F, f - p.data))
	log.s = sample(log, n, replace=TRUE)
	return(M[log.s,])
}

write.predictions = function(data, u.prb, p.prb, prior, filename="uORF_function_unlabeled_predictions_sorted_") {
	prob.ratios = n.Bayes.score(data, u.prb, p.prb, prior)
	uORFs.ratio = data.frame(prob.ratios)
	rownames(uORFs.ratio) <- rownames(data)
	uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame

	b = uORFs.ratio[order(-uORFs.ratio[,1]),]
	sorted.scores = data.frame(b[,1])
	rownames(sorted.scores) <- rownames(b)
	
	output = paste(filename, sep="")
	write.table(sorted.scores, file=output, quote=FALSE, sep="\t", col.names=FALSE)
	return(output)
}

hist.ratio = function(pred, xmin, xmax, breaks=100) {
	hist(pred[pred > xmin & pred < xmax], breaks=breaks)
	
}


#-------------------------------------

print("Loading discretization library...")

library(discretization, lib.loc="~/R/x86_64-unknown-linux-gnu-library/2.11/")

disc_file = "discretized_ATG_uORF_table_2015-03-21"

## NEW MODIFICATION, for the partitioned MDLP protocol.

setwd("~/Rfiles")

all_columns <- read.table("discretized_uORF_table_all_2", header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
all_columns <- all_columns[,c(1), drop=FALSE]

for (i in 3:89) {
column_current <- read.table(paste("discretized_uORF_table_all_", i, sep=""), header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
column_current <- column_current[,c(1), drop=FALSE]
all_columns <- cbind(all_columns, column_current)
print(paste("Binding column #", i, "...", sep=""))
}

final_column <- read.table(paste("discretized_uORF_table_all_", i, sep=""), header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
final_column <- final_column[,c(2), drop=FALSE]

print("Completing input reading...")

disc.data <- cbind(all_columns, final_column)

write.table(disc.data, file=paste("./discretized_data_all_1.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)


print("Reading discretized file...")

# disc.data = read.table(disc_file, header=TRUE,sep="\t",row.names=1,comment.char="", quote="")

# Imaging the table ------------------------------------

hist.trial = function(pred, xmin, xmax, breaks=30, main="No main provided", xlab="value", col="lightgreen") {
	hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
	
}

#The names of the features before removal of homogenous features are written in ./featuresbefore.txt

allfeaturenamesbefore = colnames(disc.data)
write.table(allfeaturenamesbefore, file=paste("./featuresbefore.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#The data from each of the discretized features (both unl and pos combined) is plotted next, on an individual basis.

xlim = c(0, 10000)
xlab = paste("trial", sep="")
tracker = 1

for (aname in allfeaturenamesbefore) {

#samplename = paste("./plots/", allfeaturenames[k], ".pdf", sep="")
pdf(paste("./plots/", aname, ".pdf", sep=""))
hist.trial(disc.data[,tracker], xlim[1], xlim[2], breaks=30, main=aname, xlab, col="lightblue")
dev.off()
tracker <- tracker + 1

}

# Imaging the table-------------------------------------

print("Removing homogenous features...")

disc.data = screen(disc.data)			# Throw out homogenous features (all 1's or all 2's)

write.table(disc.data, file=paste("./discretized_data_all_1_screened.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

#The names of each feature, retained following discretization, are written to the file ./featuresafter.txt

allfeaturenamesafter = colnames(disc.data)
write.table(allfeaturenamesafter, file=paste("./featuresafter.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#the names of the rejected features (homogenous features) are written to the file ./rejected.txt

rejected = setdiff(allfeaturenamesbefore, allfeaturenamesafter)
write.table(rejected, file=paste("./rejected.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# print("Homogenous features rejected.")

#the dimensions of the discretized data are measured

d.M = dim(disc.data)[1]
d.N = dim(disc.data)[2]

#the positive and unlabeled data sets are separated based on the class column, the class column is removed.

all_pos.disc = (disc.data[disc.data[,d.N] == 1,])[, -d.N]	# shaves off the class column
all_unl.disc = (disc.data[disc.data[,d.N] == 0,])[, -d.N]	# shaves off the class column
all_neu.disc = (disc.data[disc.data[,d.N] == 2,])[, -d.N]	# shaves off the class column

# #retaining a number of the examples from the positive set, to test the quality of the algorithm:

# d.p.M = dim(all_pos.disc)[1] #number of rows in the positive example.
# d.p.N = dim(all_pos.disc)[2] #number of columns in the positive example.
# tenpercentsamplesize = trunc(d.p.M/10)
# tenpercentremainder = d.p.M %% 10

# #In the following section of code, looping over the data set is accomplished, to make multiple retained sets, and multiple products.

# #First we randomly sample the positive data:

# random_max = length(all_pos.disc[,1])
# nonrandom_vector = seq(1,random_max,1)
# random_vector = sample(nonrandom_vector)
# all_pos.disc.mixed = all_pos.disc[random_vector,]

# retained = list()
# for (i in 0:9) {
# retained[[(i+1)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):((i+1)*tenpercentsamplesize),1:d.p.N]
# }

# if (tenpercentremainder == 0) {
# positivestail = list()
# positivestail[[10]] = list()
# for (i in 1:9) {
# positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
# }
# } else {
# positivestail = list()
# for (i in 1:10) {
# positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
# }
# }

# positiveshead = list()
# positiveshead[[(1)]] = list()
# for (i in 2:10) {
# positiveshead[[(i)]] = all_pos.disc.mixed[1:((i-1)*tenpercentsamplesize),1:d.p.N]
# }

# positives = list()
# for (i in 1:10) {
# positives[[(i)]] = rbind(positiveshead[[i]],positivestail[[i]])
# }

# #now add the positive test group, back to the unlabeled set, for later retrieval

# bound = list()
# for (i in 1:10) {
# bound[[i]] = rbind(retained[[i]],all_unl.disc)
# }

# #-----------------------------------------------
# # The following will be looped:

# #The dimensions of both the unlabeled, and positive, matrices of discretized values are measured.

# for (i in 1:10) {

# setwd("~/Rfiles/")

# unl.looping.disc = bound[[i]]
# pos.looping.disc = positives[[i]]
# ret.looping.disc = retained[[i]]

# u.M = dim(unl.looping.disc)[1]
# print(u.M)
# u.N = dim(unl.looping.disc)[2]
# print(u.N)
# p.M = dim(pos.looping.disc)[1]
# print(p.M)
# p.N = dim(pos.looping.disc)[2]
# print(p.N)
# print(d.N)

# #The number of variables is measured.

# max.len = 1
# for (k in 1:d.N) {
# 	len = length(table(disc.data[,k]))
# 	if (len > max.len) max.len <- len
# }

# print(max.len)

# #probability matrices are constructed (empty at this stage)

# u.prb = matrix(0, nrow=u.N, ncol=max.len)
# p.prb = matrix(0, nrow=u.N, ncol=max.len)
# dimnames(u.prb)[1] = dimnames(unl.looping.disc)[2]
# dimnames(p.prb)[1] = dimnames(pos.looping.disc)[2]

# #numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.

# for (k in 1:u.N) 
# 	for (j in 1:max.len) {
# 		u.prb[k, j] = mean(unl.looping.disc[,k] == j)
# 		p.prb[k, j] = mean(pos.looping.disc[,k] == j)		
# 	}

# prior = 0.61

# setwd("~/Sandbox/")

# #write.predictions

# unl.pred.file = write.predictions(unl.looping.disc, u.prb, p.prb, prior, filename=paste("uORF_function_unlabeled_predictions_sorted_discretized_looping", i, sep=''))
# pos.pred.file = write.predictions(pos.looping.disc, u.prb, p.prb, prior, filename=paste("uORF_function_positive_predictions_sorted_discretized_looping", i, sep=''))
# ret.pred.file = write.predictions(ret.looping.disc, u.prb, p.prb, prior, filename=paste("uORF_function_retained_predictions_sorted_discretized_looping", i, sep=''))

# }

# #Concatenate the files together:

# setwd("~/Sandbox/")

# cat.unl = list()
# cat.pos = list()
# cat.ret = list()

# for (i in 1:10) {
# catfile.unl = paste("uORF_function_unlabeled_predictions_sorted_discretized_looping", i, sep='')
# cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
# cat.unl = rbind(cat.unl,cat.unl.add)
# catfile.pos = paste("uORF_function_positive_predictions_sorted_discretized_looping", i, sep='')
# cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
# cat.pos = rbind(cat.pos,cat.pos.add)
# catfile.ret = paste("uORF_function_retained_predictions_sorted_discretized_looping", i, sep='')
# cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
# cat.ret = rbind(cat.ret,cat.ret.add)
# }

# hist.looping = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabeled probability ratio\n prior positive distribution is 0.61", col="lightgreen") {
# 	hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)	
# }

# #-------------------------------------

# svg('ret_pos_unl_distro.svg')
# par(mfrow = c(3,1))
# # xlim = c(1.00, 1000)
# xlim = c(-100, 100)
# hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
# hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabeled uORFs", col="lightgreen")
# hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
# dev.off()

# # Some statistics, resulting from the looping algorithm:

# quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
# quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
# quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))

# setwd("~/Rfiles/")

#-------------------------------------------------

#student's t-test the variables remaining, and plot ttest result

#A matrix is constructed, allowing for measurement of the t-test value for each variable.

ttestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
dimnames(ttestmat)[1] = dimnames(all_pos.disc)[2]
kstestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
dimnames(kstestmat)[1] = dimnames(all_pos.disc)[2]
tracker=1
for (name in allfeaturenamesafter[1:(length(allfeaturenamesafter)-1)]) {
ttest = t.test(all_pos.disc[,tracker],all_unl.disc[,tracker])
kstest = ks.test(all_pos.disc[,tracker],all_unl.disc[,tracker])
ttestmat[tracker,1] <- ttest$statistic
ttestmat[tracker,2] <- ttest$parameter
ttestmat[tracker,3] <- ttest$p.value
ttestmat[tracker,4] <- ttest$conf.int[1]
ttestmat[tracker,5] <- ttest$conf.int[2]

kstestmat[tracker,1] <- kstest$statistic
kstestmat[tracker,2] <- kstest$p.value
kstestmat[tracker,3] <- kstest$alternative
kstestmat[tracker,4] <- kstest$method
kstestmat[tracker,5] <- kstest$data.name
tracker <- tracker + 1
}

#Order the t-test matrix according to the magnitude of the t-test value, print to ./ttestmatordered.txt

ttestmatordered=ttestmat[order(abs(ttestmat[,1])),]
kstestmatordered=kstestmat[order(abs(as.numeric(kstestmat[,1]))),]
write.table(ttestmatordered, file=paste("./ttestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(kstestmatordered, file=paste("./kstestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)

#condense the tissue types variables:
kstestmatordered_GTEXcondensed = kstestmatordered
removal_vector = c("Bone.Marrow","Liver","Pituitary","Spleen","Bladder","Skin","Stomach","Lung","Nerve","Small.Intestine","Blood.Vessel","Muscle","Adipose.Tissue","Pancreas","Salivary.Gland","Esophagus","Blood","Brain","Thyroid","Fallopian.Tube","Vagina","Kidney","Prostate","Uterus","Cervix.Uteri","Colon","Breast","Heart","Testis","Ovary","Adrenal.Gland")
forcombined_average = kstestmatordered_GTEXcondensed[rownames(kstestmatordered_GTEXcondensed) %in% removal_vector, ]

combined_average = mean(as.numeric(forcombined_average[,1]))
combined_pvalue = mean(as.numeric(forcombined_average[,2]))
combined_type = forcombined_average[1,3]
combined_test = forcombined_average[1,4]
combined_data = forcombined_average[1,5]

kstestmatordered_GTEXcondensed = kstestmatordered_GTEXcondensed[!rownames(kstestmatordered_GTEXcondensed) %in% removal_vector, ]
kstestmatordered_GTEXcondensed = rbind(kstestmatordered_GTEXcondensed, "GTEX_combined" = c(combined_average, combined_pvalue, combined_type,combined_test, combined_data))
kstestmatordered_GTEXcondensed=kstestmatordered_GTEXcondensed[order(abs(as.numeric(kstestmatordered_GTEXcondensed[,1]))),]

write.table(kstestmatordered_GTEXcondensed, file=paste("./kstestmatordered_GTEXcondensed.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
#All variables t-test plot
pdf(paste("./plots/ttest/", "ttestordered", ".pdf", sep=""))
barplot(abs(ttestmatordered[,1]), main="t-test", horiz=TRUE, names.arg=c(1:dim(ttestmatordered)[1]), las=1)
dev.off()
#All variable ks-test plot
pdf(paste("./plots/kstest/", "kstestordered_original", ".pdf", sep=""))
barplot(abs(as.numeric(kstestmatordered[,1])), main="ks-test", horiz=TRUE, names.arg=c(1:dim(kstestmatordered)[1]), las=1)
dev.off()
#Top 10 variables ks-test plot
pdf(paste("./plots/kstest/", "kstestordered_top10", ".pdf", sep=""))
par(mar=c(5.1, 13 ,4.1 ,2.1))
barplot(abs(as.numeric(kstestmatordered_GTEXcondensed[1:10,1])), main="ks-test", horiz=TRUE, names.arg=tail(dimnames(kstestmatordered_GTEXcondensed)[[1]], 10), las=1)
dev.off()

#Legend for the t-test table is constructed, so that the numbers on the table can be associated with variables.

ttestlegend = ttestmatordered[,1:2]
ttestlegend[,2] = ttestmatordered[,1]
ttestlegend[,1] = 1:dim(ttestmatordered)[1]
write.table(ttestlegend, file=paste("./ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)

#---------------------------------------------------

#The dimensions of both the unlabeled, and positive, matrices of discretized values are measured.

u.M = dim(all_unl.disc)[1]
print(u.M)
u.N = dim(all_unl.disc)[2]
print(u.N)
p.M = dim(all_pos.disc)[1]
print(p.M)
p.N = dim(all_pos.disc)[2]
print(p.N)
print(d.N)

#The number of variables is measured.

max.len = 1
for (k in 1:d.N) {
	len = length(table(disc.data[,k]))
	if (len > max.len) max.len <- len
}

print(max.len)

#probability matrices are constructed (empty at this stage)

u.prb = matrix(0, nrow=u.N, ncol=max.len)
p.prb = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.prb)[1] = dimnames(all_unl.disc)[2]
dimnames(p.prb)[1] = dimnames(all_pos.disc)[2]

#numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.

for (k in 1:u.N) 
	for (j in 1:max.len) {
		u.prb[k, j] = mean(all_unl.disc[,k] == j)
		p.prb[k, j] = mean(all_pos.disc[,k] == j)		
	}

#A table with the number of instances of each value of each variable is constructed (both positive and unlabeled)

u.count = matrix(0, nrow=u.N, ncol=max.len)
p.count = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.count)[1] = dimnames(all_unl.disc)[2]
dimnames(p.count)[1] = dimnames(all_pos.disc)[2]

#A table with the frequency (percentwise) of each value of for each variable is constructed (both positive and unlabeled)

u.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
p.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.cumulfreq)[1] = dimnames(all_unl.disc)[2]
dimnames(p.cumulfreq)[1] = dimnames(all_pos.disc)[2]

for (k in 1:u.N)
	{
		varsp=table(all_pos.disc[,k])
		print(varsp)
		varsu=table(all_unl.disc[,k])
		print(varsu)
		varsall=table(c(all_unl.disc[,k],all_pos.disc[,k]))
		print(varsall)
		u.count[k,1:length(varsu)] = table(all_unl.disc[,k])
		p.count[k,1:length(varsp)] = table(all_pos.disc[,k])
		vars=table(all_unl.disc[,k])/length(all_unl.disc[,k])
		#print(vars)
		u.cumulfreq[k,1:length(varsu)] = (table(all_unl.disc[,k])/length(all_unl.disc[,k]))
		p.cumulfreq[k,1:length(varsp)] = (table(all_pos.disc[,k])/length(all_pos.disc[,k]))
	}

print(p.count)

print(p.cumulfreq)
print(u.cumulfreq)

print(p.cumulfreq[1,1:3])
print(u.cumulfreq[1,1:3])

#Plot with differential scores between categories ALL DIFFERENT PLOTS

tracker = 1

anames = rownames(p.cumulfreq)

for (k in 1:dim(p.cumulfreq)[1]) {

cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])

print(cumulfreqmatvect)
colnumber = length(cumulfreqmatvect)/2

trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)

trialmatnormalized = trialmat
for (index in 1:dim(trialmat)[2]) {
trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
}

pdf(paste("./plots/retained/", anames[k], ".pdf", sep=""))
barplot(trialmatnormalized, main=anames[k], xlab="value (L to R, 1 to 2 (or 3))", col=c("darkblue","red"), legend = c("positive", "unlabeled"))
dev.off()
tracker <- tracker + 1

}

# #ALL ON THE SAME PLOT

pdf(paste("./plots/retained/", "ALLPLOT", ".pdf", sep=""))
par(mfrow=c(2,2))
for (k in 1:dim(p.cumulfreq)[1]) {

cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
colnumber = length(cumulfreqmatvect)/2
trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)

trialmatnormalized = trialmat
for (index in 1:dim(trialmat)[2]) {
trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
}

barplot(trialmatnormalized, main=anames[k], xlab="value", col=c("darkblue","red"), legend = c("positive", "unlabeled"))
tracker <- tracker + 1

}
dev.off()

# 0.53 was a good prior in the cross validation
# 0.61 is the prior I chose to compare normal and corrected discretized predictions
# Ahh, it turns out that in my corrected discretized predictions the CDS.start.to.uORF.end variable has all values for positive as 1, thus giving 0 probability for any unlabeled uORF that doesn't meet that criteria. While this may be cause for concern and motivate using one of the Naive Bayes modifications to keep any probability from being 0, I think ti is just fine to leave out these examples. However, I will get an infinity error for dividing by this 0 in outputing the result

setwd("~/Sandbox/")

prior = 0.61

#write.predictions

unl.pred.file = write.predictions(all_unl.disc, u.prb, p.prb, prior, filename=paste("uORF_function_unlabeled_predictions_sorted_discretized_", date, sep=''))
pos.pred.file = write.predictions(all_pos.disc, u.prb, p.prb, prior, filename=paste("uORF_function_positive_predictions_sorted_discretized_", date, sep=''))
neu.pred.file = write.predictions(all_neu.disc, u.prb, p.prb, prior, filename=paste("uORF_function_neutral_predictions_sorted_discretized_", date, sep=''))

# ret.pred.file = write.predictions(retpos.disc, u.prb, p.prb, prior, filename=paste("uORF_function_retained_predictions_sorted_discretized_", date, sep=''))

unl.pred.file = paste("uORF_function_unlabeled_predictions_sorted_discretized_", date, sep='')
pos.pred.file = paste("uORF_function_positive_predictions_sorted_discretized_", date, sep='')
neu.pred.file = paste("uORF_function_neutral_predictions_sorted_discretized_", date, sep='')

# ret.pred.file = paste("uORF_function_retained_predictions_sorted_discretized_", date, sep='')
   
#-------------------------------------

# Function Definitions 

hist.ratio = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabeled probability ratio\n prior positive distribution is 0.61", col="lightgreen") {
	hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
	
}

#-------------------------------------
setwd("~/Sandbox/")

unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="")
pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="")
neu.pred = read.table(neu.pred.file, header=FALSE,sep="\t",comment.char="", quote="")

unl.pred.samp = unl.pred[sample(1:length(as.matrix(unl.pred[,1])),10000),]

library("ROCR")
target_pred = rbind(pos.pred,unl.pred.samp)
pos.neu.pred = rbind(pos.pred,neu.pred)
class.pos <- matrix(sample(1, (ncol(pos.pred)*nrow(pos.pred)), replace=T), ncol=ncols)
class.unl <- matrix(sample(0, (ncol(unl.pred.samp)*nrow(unl.pred.samp)), replace=T), ncol=ncols)
target_class <-rbind(class.pos,class.unl)
pred <- prediction(target_pred, target_class)
perf <- performance(pred,"tpr","fpr")
pred_comb <- prediction(target_pred_comb, target_class_comb)
perf_comb <- performance(pred,"tpr","fpr")

require(ggplot2)

#This ROC curve only counts the positive data (excludes neutral from the analysis).

svg("./ROC_prime.svg")
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
plot(perf,col="black",lty=3, lwd=3)
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc<-round(auc, digits = 2)
auct <- paste(c("(AUC)  = "),auc,sep="")
legend(0.3,0.6,c(auct,"\n"),border="white",cex=1.7,box.col = "white")
dev.off()

print(pos.pred[1,2])
print(unl.pred[1,2])
print(neu.pred[1,2])

bench_hist = hist(unl.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabeled uORFs", col="lightblue")

svg('rplot.svg')
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist(unl.pred[,2], xlim[1], xlim[2], breaks=bench_hist$breaks, main="Class probability ratio distribution for unlabeled uORFs", col="lightblue")
hist(pos.pred[,2], xlim[1], xlim[2], breaks=bench_hist$breaks, main="Class probability ratio distribution for positive uORFs")
hist(neu.pred[,2], xlim[1], xlim[2], breaks=bench_hist$breaks, main="Class probability ratio distribution for neutral uORFs", col="lightblue")
dev.off()

svg('rplot_neu.svg')
par(mfrow = c(1,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
dev.off()

svg('rplot2.svg')
par(mfrow = c(2,1))
# xlim = c(1.00, 1000)
xlim = c(0, 100)
hist(unl.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabeled uORFs\n predicted as positives", col="lightblue")
hist(pos.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs\n predicted as positives")
dev.off()

svg('rplot2_neu.svg')
par(mfrow = c(1,1))
# xlim = c(1.00, 1000)
xlim = c(0, 100)
hist(neu.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for neutral uORFs\n predicted as positives", col="lightblue")
dev.off()

high.thresh = quantile(unl.pred[,2], probs = c(0.93))
high.preds = unl.pred[unl.pred[,2] >= high.thresh,]

write.table(high.preds, file="top7%_predicted_functional_uORFs.txt", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

print(quantile(pos.pred[pos.pred[,2] > 1.00, 2], probs = c(0.45)))
print(quantile(unl.pred[unl.pred[,2] > 1.00, 2], probs = c(0.45)))
summary(pos.pred[pos.pred[,2] > 1.00, 2])

summary(unl.pred[unl.pred[,2] > 1.00, 2])

#------------------------------------------

# Distribution plots and Box and Whisker plots to evaluate positive predictions

hist.ratio = function(pred, xmin, xmax, breaks=100, label="unlabeled", col='lightblue', prior=0.61) {
	hist(pred[pred > xmin & pred < xmax], xlim=c(xmin, xmax), col=col, breaks=breaks, main=paste("Class probability ratio distribution\n for", label, "uORFs\n predicted as positives"), xlab=paste("Positive/unlabeled probability ratio\n with", prior, "prior distribution"))
	
}

get.pos = function(data) {
	return(data[data[,2] > 1 & !is.na(data[,2]),])
}

# Discretized results with a 0.61 prior.

setwd("~/Sandbox/")

# disc.file = "uORF_function_unlabeled_predictions_sorted_discretized_2014-09-23"

# disc.pred = read.table(disc.file, header=FALSE,sep="\t",comment.char="", quote="")

# disc.pos.file = "uORF_function_positive_predictions_sorted_discretized_2014-09-23"

# disc.pos = read.table(disc.pos.file, header=FALSE,sep="\t",comment.char="", quote="")

d.unl = get.pos(unl.pred)
d.pos = get.pos(pos.pred)

svg('rplot3.svg')
boxplot(list(d.unl[,2], d.pos[,2]) , outline=FALSE, boxwex=0.3, names=c("Unlabeled", "Positive"), main="Class probability ratio quartiles\n for uORFs\n predicted as positives")
# boxplot(d.pos[,2], outline=FALSE, add=TRUE)
dev.off()
svg('rplot4.svg')
boxplot(list(unl.pred[,2], pos.pred[,2]) , outline=FALSE, boxwex=0.3, names=c("Unlabeled", "Positive"), main="Class probability ratio quartiles\n for uORFs")
# boxplot(d.pos[,2], outline=FALSE, add=TRUE)
dev.off()

comb.pred = rbind(unl.pred,pos.pred,neu.pred)
comb.pred.finite = (comb.pred[is.finite(comb.pred[,2]),])

sampling_number = min(c(length(neu.pred[,2]), length(pos.pred[,2])))

pos.dist.funct = ecdf(sample(pos.pred[,2], sampling_number, replace = FALSE))
unl.dist.funct = ecdf(sample(unl.pred[,2], sampling_number, replace = FALSE))
neu.dist.funct = ecdf(sample(neu.pred[,2], sampling_number, replace = FALSE))
cumul.dist.funct = ecdf(sample(comb.pred[,2], sampling_number, replace = FALSE))

svg('pos_unl_combined_CDF.svg')
par(mfrow = c(3,1))
plot(unl.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='green', pch="", lwd=10)
plot(pos.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='red', pch="", lwd=10)
plot(neu.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='blue', pch="", lwd=10)

dev.off()

xandytest = pos.dist.funct

# plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF

xandytest = as.list(environment(pos.dist.funct))

x = xandytest[[2]]
y = xandytest[[3]]

print(x)
print(y)
print(length(x))
print(length(y))

svg('pos_CDF.svg')
spl = smooth.spline(x, y, df = 10)

pred = predict(spl, deriv = 1)
normpredy = pred$y / max(pred$y)

pred1 = predict(spl, deriv = 2)
normpred1y = pred1$y / abs(min(pred1$y))

plot (x, y, ylim=c(-1, 1))
lines(spl, col=1)
lines(pred$x, normpredy, col=2)
lines(pred1$x, normpred1y, col=3)

legend(-20,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))

dev.off()


# plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF

xandytest = as.list(environment(cumul.dist.funct))

x = xandytest[[2]]
y = xandytest[[3]]

print(x)
print(y)
print(length(x))
print(length(y))

svg('unl_CDF.svg')
spl = smooth.spline(x, y, df = 10)

pred = predict(spl, deriv = 1)
normpredy = pred$y / max(pred$y)

pred1 = predict(spl, deriv = 2)
normpred1y = pred1$y / max(pred1$y)

plot (x, y, ylim=c(-1, 1))
lines(spl, col=1)
lines(pred$x, normpredy, col=2)
lines(pred1$x, normpred1y, col=3)

legend(10,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))

dev.off()


print(xandytest[3])
print(length(comb.pred[,2]))
x = xandytest[[2]]
y = xandytest[[3]]
xandytest[3]
print(length(xandytest))


svg('rplot7.svg')
ycs.prime = diff(xandytest)/diff(comb.pred[,2])
pred.prime = predict(spl, deriv=1)

# plot(ycs.prime)
# lines(pred.prime$y, col=2)


# dev.off()

# svg('rplot6.svg')
# cumul.dist.funct.deriv = predict(c(comb.pred[,2],xandytest), cumul.dist.funct[,1], 1)
# plot(c(comb.pred[,2],xandytest), xlim=c(-30, 30), main="Cumulative Distribution Derivative Function for all uORFs")
# dev.off()
