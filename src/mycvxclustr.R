rm(list=ls())
setwd("/home/duma/TCGA/xpca")

#source("create_adjacency.R")
#source("find_clusters.R")

library(igraph)
library(Matrix)
library(compiler)
library(cvxclustr)

wd = "bc"

#args = commandArgs(trailingOnly = TRUE);
#print(args);
#k = as.numeric(args[1])

## Read in data matrix (svd/mds/nmf/xpca)
X = read.csv(paste("../data/", wd, "/U.csv", sep=""), row.names = 1) 

##DEBUG
print("str(X)")
print(str(X))
rows = rownames(X)
cols = colnames(X)

X = as.matrix(X)
X = X/norm(as.matrix(X),'f')
#X = t(scale(X,center=TRUE,scale=FALSE))
X = t(X)
print("dim of matrix X")
print(dim(X))
m = nrow(X)
n = ncol(X)

## Pick some weights and a sequence of regularization parameters.
k = 10
phi = 0.5
w = kernel_weights(X, phi)
w = knn_weights(w, k, n)
print("dim and type of w")
#print(w[w < 1e-10])
#print(str(w))

print("weights graph")
A = weights_graph(w, n)
print(find_clusters(A))

print("dim and type of wh")
wh = which(w > 0)
print(length(wh))
print(str(wh))

nk = n*(n-1)/2
print("nk")
print(nk)

#gamma = seq(0.0, 43, length.out = 100)
#gamma = seq(0, 101, by = 1)
gamma = 10^seq(-4,0,length.out=1e2)


# Perform clustering
nu = AMA_step_size(w,n)
print("nu")
print(nu)
print("1/n")
print(1/n)
sol = cvxclust(X, w, gamma, nu, method = "ama") # nu = 1.999/n

print("summary(sol)")
print(summary(sol))

## vec2tri
vec2tri <- cmpfun(function(k,p) {
i <- ceiling(0.5*(2*p-1 - sqrt((2*p-1)^2 - 8*k)))
j <- k - p*(i-1) + i*(i-1)/2 + i
return(as.matrix(cbind(i,j)))
})

nk = n*(n-1)/2
ix = vec2tri(which(w > 0), n)
print ("dim and type of ix")
print(dim(ix))
print(str(ix))

nGamma = sol$nGamma
for (j in 1:nGamma) {
	## Construct adjacency matrix
	A = create_adjacency(sol$V[[j]], ix, n, method = "ama")
	#print ("dim and type of A")
	#print(dim(A))
	#print(str(A))
	fclusters = find_clusters(A)
	#print("fclusters type")
	#print(str(fclusters))
	print("numbers of clusters")
	print(length(fclusters$size))
	#print("head of clusters")
	#print(head(fclusters$size))
}

############################Plot the cluster path#############################
#library(ggplot2)

#svdX = svd(X)

#pc = svdX$u[,1:2,drop=FALSE]
#print("type of pc")
#print(str(pc))
#pc.df = as.data.frame(t(pc)%*%X)
#print("type of pc.df")
#print(str(pc.df))

#nGamma = sol$nGamma

#df.paths = data.frame(x=c(),y=c(), group=c())

#for (j in 1:nGamma) {
#	pcs = t(pc)%*%sol$U[[j]]
#	x = pcs[1,]
#	y = pcs[2,]
#	df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
#	df.paths = rbind(df.paths,df)
#}
#print("dim of df.paths")
#print(dim(df.paths))

#X_data = as.data.frame(t(X)%*%pc)
#print("dim and type of X_data")
#print(dim(X_data))
#print(str(X_data))
#colnames(X_data) = c("x","y")
#X_data$Name = rows

#data_plot = ggplot(data=df.paths,aes(x=x,y=y))
#data_plot = data_plot + geom_path(aes(group=group),colour="grey30",alpha=0.5)
#data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name), position=position_jitter(h=0.125,w=0.125))
#data_plot = data_plot + geom_point(data=X_data, aes(x=x,y=y), size=2)
#data_plot = data_plot + xlab("Principal Component 1") + ylab("Principal Component 2")
#data_plot + theme_bw()
#print(data_plot)

#print to file
#outFile = "../bc_mrna_figs/svd/cvx_initU.eps"
#dev.print(device=postscript, outFile, onefile=FALSE, horizontal=FALSE)
