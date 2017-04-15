rm(list=ls())
setwd("/home/duma/TCGA/firehose/bc_mut/gdac.broadinstitute.org_BRCA.Mutation_Packager_Calls.Level_3.2014101700.0.0/")

filenames = list.files(".", pattern="*.txt", full.names=TRUE)
print("head(filenames)") 
print(head(filenames)) 
print("length(filenames)")
print(length(filenames))
#ldf = lapply(filenames, function(l) read.table(l, header = TRUE, sep = "\t"))
#names(ldf) = filenames
#print("str(ldf)")
#print(str(ldf))

for (i in 1:length(filenames)) {
	print(filenames[[i]])
	df = read.table(filenames[[i]], header = FALSE, sep = "\t")
	print(dim(df))
}

#for (i in 1:length(ldf)) {
#	print(dim(ldf[[i]]))
#}
