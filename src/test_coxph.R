setwd("/home/duma/TCGA/xpca")

library(survival)
rm(list=ls())
hmohiv = read.table("http://www.ats.ucla.edu/stat/r/examples/asa/hmohiv.csv", sep=",", header = TRUE) 
print(str(hmohiv))
attach(hmohiv)
age.coxph = coxph( Surv(time,censor)~age, method="breslow")
print(summary(age.coxph))
