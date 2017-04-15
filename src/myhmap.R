rm(list=ls())

setwd("/home/duma/TCGA/xpca")

# read in data matrix
U = scan(file="Y.csv", what=character());
U.ls = lapply(U, function(l) unlist(strsplit(l, ",")))
cname = sapply(U.ls[-1], function(l) l[1])
vname = U.ls[[1]]
print("length of U.ls, cname and vname");
print(length(U.ls))
print(length(cname))
print(length(vname))

U.mat = sapply(2:length(U.ls), function(i) as.numeric(U.ls[[i]][-1]))
#U.mat = t(U.mat)
print("dim and type of U.mat")
print(dim(U.mat))
print(str(U.mat))
colnames(U.mat) = cname
rownames(U.mat) = vname

# read in patient labels
types = readLines("../data/bc_mRNA_types_types.txt", n = -1);
#unlink("../data/bc_mRNA_types_types.txt");
print("types type and size");
print(str(types));

#cc = rainbow(ncol(U.mat))
#patient lables (cancer subtypes)
color.map = function(type) { 
	if (type=="LuminalA") { 
		"blue" 
	} else if (type=="LuminalB") {
		"cyan" 
	} else if (type=="HER2") {
		"red" 
	} else if (type=="Basal") {
		"green"
	} else if (type=="Normal") {
		"black"
	}
}
#print("color.map type and size")
#print(str(color.map))

patientcolors = unlist(lapply(types, color.map))
print("patientcolors type and size")
print(str(patientcolors))

#create heatmap
#heatmap(U.mat, Rowv = NA, col = topo.colors(100), ColSideColors = patientcolors, cexCol = 0.5);

#create heatmap2
library("gplots")
f = function(d, method = "average", members = NULL) {
	hclust(d, method=method, members)
}

#postscript(file="initU.eps")
heatmap.2(U.mat, dendrogram = "column", scale="row", hclust = f, col=redgreen(75), ColSideColors=patientcolors, key=TRUE, density.info="none", trace="none")
#dev.off()
dev.print(device=postscript, "Y.eps", onefile=FALSE, horizontal=FALSE)

