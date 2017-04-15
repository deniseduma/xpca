#!/usr/bin/Rscript

rm(list=ls())
setwd("/home/duma/TCGA/xpca")

library(pamr)
#library(survival) 

args=commandArgs(TRUE)
print("str(args)")
print(str(args))

data = args[1] #"met"
inputdir = args[2] #"bc_met"
method = args[3] #"svd"
k = args[4] #2

classifier = "K"
if (method == "svd") {
	title = paste("initU_", classifier, "_k", k,  sep = "") 
} else if (method == "mds") {
	title = paste("Y_", classifier, "_k", k, sep = "") 
} else if (method == "nmf") {
	title = paste("Z_", classifier, "_k", k, sep = "") 
} else if (method == "xpca") {
	title = paste("U_", classifier, "_k", k, sep = "") 
}

## Read CSV
## BC/OV/UCEC/LUAD/GBM
#inFile = "../data/bc_mut_cna_clusters_U.csv"
inFile = paste("../data/", inputdir, "/", title, ".csv", sep = "")
print(paste("inFile", inFile, sep=" "))
U = read.table(inFile, sep=",", header = TRUE) 
print("str(U)")
print(str(U))

#pamr.plotsurvival(as.factor(U$Cluster), U$OS_Time, U$OS_Event)
cluster.coxph = coxph( Surv(U$OS_Time,U$OS_Event)~as.factor(U$Cluster))
print(summary(cluster.coxph))

## Save eps to file
# BC / OV / LUAD / GBM / UCEC
#outFile = "../bc_mrna_figs/xpca/survival_U_K.eps"
#outFile = paste("../", inputdir, "_figs/", method, "/", title, ".eps", sep="")
#dev.print(device=postscript, outFile, onefile=FALSE, horizontal=FALSE)


