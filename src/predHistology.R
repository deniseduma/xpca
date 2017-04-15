#!/usr/bin/Rscript

rm(list=ls())
setwd("/home/duma/TCGA/xpca/")

#############################EXAMPLE########################
#set.seed(123)
#N        <- 50
#smokes   <- factor(sample(c("no", "yes"), N, replace=TRUE))
#siblings <- factor(round(abs(rnorm(N, 1, 0.5))))
#cTab     <- table(smokes, siblings)
#addmargins(cTab)

##	siblings
##smokes  0  1  2 Sum
##     no   5 16  4  25
##     yes  3 19  3  25
##     Sum  8 35  7  50

##chisq.test(cTab)

#chisq.test(cTab)

##Pearson's Chi-squared test
##data:  cTab 
##X-squared = 0.9, df = 2, p-value = 0.6376

args=commandArgs(TRUE)
print("str(args)")
print(str(args))

data = args[1] #"met"
inputdir = args[2] #"ucec_met"
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

inFile = paste("../data/", inputdir, "/", title, ".csv", sep = "")
U = read.table(inFile, sep=",", header = TRUE) 
print("str(U)")
print(str(U))

V = table(factor(U$Histology), factor(U$Cluster))
addmargins(V)

print("type of V")
print(str(V))
print("V")
print(V)

test = chisq.test(V)

print(test)

