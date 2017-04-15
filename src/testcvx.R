rm(list=ls())

setwd("/home/duma/TCGA/xpca")

library(cvxclustr)

## Clusterpaths for Mammal Dentition
data(mammals)
print("type of mammals")
print(mammals)
X = as.matrix(mammals[,-1])
print("type of orig X")
print(dim(X))
X = t(scale(X,center=TRUE,scale=FALSE))
p = ncol(X)
print("type of trans X")
print(dim(X))

## Pick some weights and a sequence of regularization parameters.
k = 5
phi = 0.5
w = kernel_weights(X,phi)
w = knn_weights(w,k,p)
gamma = seq(0.0,43, length.out=100)

## Perform clustering
sol = cvxclust(X, w, gamma, nu = 1.999/p)
print("type of sol")
print(length(sol))

### Plot the cluster path
library(ggplot2)

svdX = svd(X)

pc = svdX$u[,1:2,drop=FALSE]
print("dim and type of pc")
print(dim(pc))
print(str(pc))
pc.df = as.data.frame(t(pc)%*%X)
print("dim and type of pc.df")
print(dim(pc.df))
#print(str(pc.df))

nGamma = sol$nGamma

df.paths = data.frame(x=c(),y=c(), group=c())

for (j in 1:nGamma) {
	pcs = t(pc)%*%sol$U[[j]]
	x = pcs[1,]
	y = pcs[2,]
	df = data.frame(x=pcs[1,], y=pcs[2,], group=1:p)
	df.paths = rbind(df.paths,df)
	#print("dim of df.paths")
	#print(dim(df.paths))
}

print("dim and type of df.paths")
print(dim(df.paths))
#print(str(df.paths))

X_data = as.data.frame(t(X)%*%pc)
colnames(X_data) = c("x","y")
X_data$Name = mammals[,1]
print("dim and type of X_data")
print(dim(X_data))
print(X_data)

data_plot = ggplot(data=df.paths, aes(x=x,y=y))
data_plot = data_plot + geom_path(aes(group=group), color="grey30", alpha=0.5)
data_plot = data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name), position=position_jitter(h=0.25,w=0.25))
data_plot = data_plot + geom_point(data=X_data, aes(x=x,y=y), size=3)
data_plot = data_plot + xlab("Principal Component 1") + ylab("Principal Component 2")
data_plot + theme_bw()

print(data_plot)
