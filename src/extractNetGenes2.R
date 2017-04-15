rm(list=ls())
setwd("/home/duma/TCGA/xpca")

wd = 'bc'

geneFile = "../data/STRING.txt"
inFile = paste("../data/", wd, "/", wd,"_mRNA_tum_small_firehose.txt", sep = "")
outFile = paste("../data/", wd, "/", wd,"_mRNA_tum_small_net_firehose.txt", sep = "")

U = read.table(inFile, sep=" ", header = TRUE, row.names = 1) 
print("dim(U)")
print(dim(U))

GENES = readLines(geneFile)
GENES = strsplit(GENES, "\t")
print("Number of interactions")
print(length(GENES))
print("Example interactions")
print(head(GENES))

V1 = sapply(GENES, function(l) l[1])
V2 = sapply(GENES, function(l) l[2])
#print("str(V1)")
#print(str(V1))
#print("str(V2)")
#print(str(V2))
VV = c(V1, V2)
#print("str(VV)")
#print(str(VV))
I = unique(VV)
print("Number of unique genes")
print(length(I))

#U = U[, which(GENES %in% colnames(U))]
U = U[, which(colnames(U) %in% I)]
print("dim(U)")
print(dim(U))

write.table(U, file = outFile, quote = FALSE)
