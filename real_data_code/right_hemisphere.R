library(R.matlab)
source('../lfmm_main.R')
library(cate)
library(dSVA)

y_right=readMat('y_right_rad.mat')
y_right=y_right[[1]]


x=read.table('x.txt')
x1=x[,1:8]
colnames(x1)=c("gender","hand","edu","age","MCI","AD","ADNIgo","ADNI2")
x1=scale(x1)
y_right=scale(y_right)

##detect number of confounders
#est.confounder.num(~ . | 1, data.frame(x1), y_right, method = "ed")
#est.confounder.num(~ .| 1, data.frame(x1), y_right, method = "bcv")

k=3
##confounder results
lfmm.result<-lfmm(y=y_right, x=x1, k=k)

##beta results
####lfmm


x2= model.matrix(~1 + x1)
lfmm_refit <- dSVA:::Get_Pvalue(y_right, x2, lfmm.result)
sum(abs(lfmm_refit$Pvalue[3,])<0.05)
sum(abs(lfmm_refit$Pvalue[4,])<0.05)
sum(abs(lfmm_refit$Pvalue[5,])<0.05)
sum(abs(lfmm_refit$Pvalue[6,])<0.05)

###dsva
dSVA.result<-dSVA(y_right, x1, ncomp=k)
#pvalue.raw<-F.mp(y_right~x1+dSVA.result$Z,which=4)$pvalue
re <- dSVA:::Get_Pvalue(y_right, x2, dSVA.result$Z)
sum(abs(re$Pvalue[3,])<0.05)
sum(abs(re$Pvalue[4,])<0.05)
sum(abs(re$Pvalue[5,])<0.05)
sum(abs(re$Pvalue[6,])<0.05)


####leap-RR

rr.result<-cate.fit(X.primary = x1, X.nuis = NULL,Y=y_right, r=k, 
                    fa.method ="ml", adj.method = "rr", 
                    psi = psi.huber, calibrate = F)

leap_refit <- dSVA:::Get_Pvalue(y_right, x2, rr.result$Z)
sum(abs(leap_refit$Pvalue[3,])<0.05)
sum(abs(leap_refit$Pvalue[4,])<0.05)
sum(abs(leap_refit$Pvalue[5,])<0.05)
sum(abs(leap_refit$Pvalue[6,])<0.05)


####unadjust
uadj <- dSVA:::Get_Pvalue(y_right, x2, NULL)
sum(abs(uadj$Pvalue[3,])<0.05)
sum(abs(uadj$Pvalue[4,])<0.05)
sum(abs(uadj$Pvalue[5,])<0.05)
sum(abs(uadj$Pvalue[6,])<0.05)


##save pvalues
p_edu=cbind(-log10(lfmm_refit$Pvalue[3,]),-log10(leap_refit$Pvalue[3,]),-log10(re$Pvalue[3,]),-log10(uadj$Pvalue[3,]))
colnames(p_edu)=c("LFMM","LEAP","DSVA","UNADJ")
write.csv(p_edu,file = "../p_edu_right.csv",row.names = F)

p_age=cbind(-log10(lfmm_refit$Pvalue[4,]),-log10(leap_refit$Pvalue[4,]),-log10(re$Pvalue[4,]),-log10(uadj$Pvalue[4,]))
colnames(p_age)=c("LFMM","LEAP","DSVA","UNADJ")
write.csv(p_age,file = "../p_age_right.csv",row.names = F)


p_ad=cbind(-log10(lfmm_refit$Pvalue[6,]),-log10(leap_refit$Pvalue[6,]),-log10(re$Pvalue[6,]),-log10(uadj$Pvalue[6,]))
colnames(p_ad)=c("LFMM","LEAP","DSVA","UNADJ")
write.csv(p_ad,file = "../p_ad_right.csv",row.names = F)


