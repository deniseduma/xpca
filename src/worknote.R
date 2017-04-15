load('unc.edu_IlluminaHiSeq_RNASeqV2.normalized.Counts.Max.saved')
colnames(mrna1.cln.max) = substr(colnames(mrna1.cln.max),1,12)

subtypes= read.table('tcga-brca-subtypes.txt',header=T,sep="\t",stringsAsFactors=F,row.names=1)

index = intersect(colnames(mrna1.cln.max),rownames(subtypes))

scolor = subtypes[index,'ER.Status']
scolor[scolor=='Positive']="red"
scolor[scolor=='Negative']="green"
scolor[is.na(scolor)]="black"
scolor[scolor=='Not Performed']="grey"
scolor[scolor=='Performed but Not Available']="pink"
scolor[scolor=='Indeterminate']="blue"

plot(log(1+as.numeric(mrna1.cln.max['ESR1',index])) , col=scolor)

rsd = apply(mrna1.cln.max,1,sd)
gindex = names(rank(rsd)[1:500])

myPCA(log(1+mrna1.cln.max[gindex,index]),label=scolor,pch=16)
############

index2 = intersect(index,rownames(subtypes[!is.na(subtypes[,'PAM50.mRNA']),]))

pcolor = subtypes[index2,'PAM50.mRNA']
pcolor[is.na(pcolor)]="white"
pcolor[pcolor=="Luminal A"]="blue"
pcolor[pcolor=="Luminal B"]="green"
pcolor[pcolor=="Basal-like"]="red"
pcolor[pcolor=="HER2-enriched"]="pink"
pcolor[pcolor=="Normal-like"]="black"

mp=myPCA(log(1+mrna1.cln.max[gindex,index2]),label=pcolor,shape=16)
mp=myPCA(mrna1.cln.max[gindex,index2],label=pcolor,shape=16)

plot(mp[,1],mp[,2],col=pcolor,pch=16,cex=1)
plot(mp[,2],mp[,3],col=pcolor,pch=16,cex=1)

write.table(mrna1.cln.max[gindex,index2],file="mRNA-seq-v2.txt",sep="\t")
write.table(pcolor,file="subtypes",sep="\t")

d = dist(mrna1.cln.max[gindex,index2])
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",main="Metric MDS",col=pcolor,pch=16)


