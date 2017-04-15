rm(list=ls())
setwd("/home/duma/TCGA/xpca")

library(NMF)

k = 5    
inputdir = "ucec_mut"

inFile = paste('../data/' , inputdir , '/', inputdir, '.txt', sep="")
print(inFile)
X = read.table(inFile, sep=" ", header = TRUE, row.names = 1) 
print("dim(X)")
print(dim(X))

#do nmf
res = nmf(X, k, nmfAlgorithm("brunet"), seed = "nndsvd", .options='v')
print("summary(res)")
print(summary(res))

# estimated target matrix
#V.hat = fitted(res)
#print("dim(V.hat)")
#print(dim(V.hat))

# Retrieve basis matrix (metagenes)
W = basis(res)
W = matrix(as.numeric(sprintf("%.4f", W)), nrow = nrow(W), ncol = ncol(W))
print("dim(W)")
print(dim(W))
# Retrieve mixture coefficient matrix (metaprofiles)
#H = coef(res)
#print("str(H)")
#print(dim(H))

# Predict row (sample) clusters
#clusters = predict(res, 'rows')
##C = matrix(as.numeric(clusters))
#print("str(clusters)")
#print(str(clusters))
#print(head(clusters))

# Print NMF decomposition to file
cols = numeric(length=k)
for (i in 1:k ) {
	cols[i] = paste("C", i, sep = "")
}
colnames(W) = seq(1:k)
rownames(W) = rownames(X)
#rownames(C) = rownamies(X)

## Print NMF decomposition to file
outFile = paste('../data/', inputdir, '/', inputdir, '_nmf_k', k, '.txt', sep="")
write.table(W, outFile, sep=" ", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Print NMF clustering 
#outFile2 = "../data/ov_mRNA_types.txt"
#write.table(C, outFile2, sep=" ", row.names = TRUE, col.names = FALSE, quote = FALSE)

