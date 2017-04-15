rm(list=ls())
setwd("/home/duma/TCGA/xpca/")

library(igraph)
library(Matrix)
library(compiler)
library(cvxclustr)

wd = "bc"
inFile = paste("../data/", wd, "/U.csv", sep = "")
U = read.csv(inFile, header = TRUE) 

X <- as.matrix(U[,-1])
X <- t(scale(X,center=TRUE,scale=FALSE))
X <- X/norm(as.matrix(X),'f')
n <- ncol(X)

k <- 3
phi <- 1e3
w <- kernel_weights(X,phi)
w <- knn_weights(w,k,n)
gamma <- 10^seq(-4,0.3, length.out=300)

## Perform clustering
nu <- AMA_step_size(w,n)
sol <- cvxclust_path_ama(X,w,gamma,nu=nu)
print("summary(sol)")
print(summary(sol))

## Plot the cluster path
library(ggplot2)

#print to file
outFile = paste("../data/", wd, "/U_Clusters.pdf", sep="")
pdf(outFile)

svdX <- svd(X)
pc <- svdX$u[,1:2,drop=FALSE]
pc.df <- as.data.frame(t(pc)%*%X)
nGamma <- sol$nGamma
df.paths <- data.frame(x=c(),y=c(), group=c())
for (j in 1:nGamma) {
  pcs <- t(pc)%*%sol$U[[j]]
  x <- pcs[1,]
  y <- pcs[2,]
  df <- data.frame(x=pcs[1,], y=pcs[2,], group=1:n)
  df.paths <- rbind(df.paths,df)
}
X_data <- as.data.frame(t(X)%*%pc)
colnames(X_data) <- c("x","y")
X_data$Name <- U[,1]
data_plot <- ggplot(data=df.paths,aes(x=x,y=y))
data_plot <- data_plot + geom_path(aes(group=group),colour='grey30',alpha=0.5)
#data_plot <- data_plot + geom_text(data=X_data,aes(x=x,y=y,label=Name),
#                                   position=position_jitter(h=0.125,w=0.125))
data_plot <- data_plot + geom_point(data=X_data,aes(x=x,y=y),size=1.5)
data_plot <- data_plot + xlab('Principal Component 1') + ylab('Principal Component 2')
data_plot + theme_bw()
plot(data_plot)

dev.off()

## Output Cluster Assignment at 235th gamma: 5 clusters
for (j in 1:nGamma) {
	A <- create_adjacency(sol$V[[j]],w,n)
	fclusters = find_clusters(A)
	print(paste("gamma", j, "numbers of clusters", length(fclusters$size), sep = " "))
	#print("head of clusters")
	#print(head(fclusters$size))
}	

