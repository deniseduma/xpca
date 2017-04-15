#!/usr/bin/Rscript

rm(list=ls())
setwd("/home/duma/TCGA/xpca")

library(NMF)
library(doParallel)
library(foreach)

args=commandArgs(TRUE)
print("str(args)")
print(str(args))

k = as.numeric(args[1])    
inputdir = args[2] #"ucec_mut"

classifier = "K"

#Read in data matrix to decompose
inFile = paste('../data/' , inputdir , '/', inputdir, '.txt', sep="")
print(inFile)
X = read.table(inFile, sep=" ", header = TRUE, row.names = 1) 
print("dim(X)")
print(dim(X))

#Read in survival data
if (inputdir == "ucec_mrna" | inputdir == "ucec_mut" | inputdir == "ucec_met") {
	suffix = "histology"
} else {
	suffix = "survival"
}
inFile2 = paste('../data/' , inputdir , '/', inputdir, "_", suffix, ".txt", sep="")
print(inFile2)
SRV = read.table(inFile2, sep=" ", header = TRUE, row.names = 1) 
print("dim(SRV)")
print(dim(SRV))

#do nmf
#res = nmf(X, k, nmfAlgorithm("brunet"), seed = "nndsvd", .options=list(verbose=TRUE, track=TRUE))
#plot(res)
res = nmf(t(X), k, nmfAlgorithm("brunet"), seed="random", nrun = 50, .pbackend=10, .options=list(parallel.required = TRUE, verbose = TRUE))
#print("summary(res)")
#print(summary(res))
## Plot a heatmap of the consensus matrix
#consensusmap(res)

# estimated target matrix
#V.hat = fitted(res)
#print("dim(V.hat)")
#print(dim(V.hat))

# Retrieve basis matrix (metagenes)
W = basis(res)
#W = matrix(as.numeric(sprintf("%.4f", W)), nrow = nrow(W), ncol = ncol(W))
print("dim(W)")
print(dim(W))
# Retrieve mixture coefficient matrix (metaprofiles)
H = coef(res)
print("dim(H)")
print(dim(H))

# Predict row (sample) clusters
clusters = predict(res, 'consensus')
print("str(clusters)")
print(str(clusters))
C = as.data.frame(clusters)
rownames(C) = rownames(X)
colnames(C) = c("Cluster")
print(paste("cophcor", cophcor(res), sep=" "))

# Print NMF decomposition to file
#cols = numeric(length=k)
#for (i in 1:k ) {
#	cols[i] = paste("C", i, sep = "")
#}
H = t(H)
colnames(H) = seq(1:k)
rownames(H) = rownames(X)

PTS = intersect(rownames(C), rownames(SRV))
C = C[PTS, ]
SRV = SRV[PTS, ]
SRVC = cbind(SRV, C)
if (inputdir == "ucec_mrna" | inputdir == "ucec_mut" | inputdir == "ucec_met") {
	colnames(SRVC) = c("Histology", "Cluster")
} else {
	colnames(SRVC) = c("OS_Event", "OS_Time", "Cluster")
}
print("dim(SRVC)")
print(dim(SRVC))

## Print NMF consensus clustering to file
outFile2 = paste("../data/", inputdir, "/Z_", classifier, "_k", k, ".csv", sep="")
write.table(SRVC, outFile2, sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Print NMF decomposition to file
outFile = paste('../data/', inputdir, '/', inputdir, '_nmf_k', k, '.txt', sep="")
write.table(H, outFile, sep=" ", row.names = TRUE, col.names = TRUE, quote = FALSE)

