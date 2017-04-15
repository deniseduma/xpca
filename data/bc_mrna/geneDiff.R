rm(list=ls())
setwd("/home/duma/TCGA/data/bc/")

l1 = readLines("bc_genes_net_sorted1.txt")
l2 = readLines("bc_genes_net_sorted2.txt")

extraGenes = setdiff(l1, l2)

writeLines(extraGenes, "bc_genes_net_extra.txt")

