rm(list=ls())
setwd("/home/duma/TCGA/xpca")

wd = 'bc'

inFile = paste("../data/", wd, "/", wd,"_mRNA_tum_small_firehose.txt", sep = "")
geneFile = paste("../data/", wd, "/", wd,"_genes_net.txt", sep = "")
outFile = paste("../data/", wd, "/", wd,"_mRNA_tum_small_net_firehose.txt", sep = "")

U = read.table(inFile, sep=" ", header = TRUE, row.names = 1) 
print("dim(U)")
print(dim(U))
GENES = readLines(geneFile)
print("length(GENES)")
print(length(GENES))
print("head(GENES)")
print(head(GENES))
print(head(colnames(U)))
print("tail(GENES)")
print(tail(GENES))
print(tail(colnames(U)))

#U = U[, which(GENES %in% colnames(U))]
U = U[, GENES]
print("dim(U)")
print(dim(U))

write.table(U, file = outFile, quote = FALSE)
