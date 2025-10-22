rm(list=ls())

source('rPCA-function-Sep16.R')
source('lfmm_main.R')

###################################### The Simulation Part######################################
n.run<-200
power<-TypeI<-AUC<-FDR<-MSE<-NUM<-NUMT<-array(NA,c(n.run,1,4,12,2))  #
Z_the<-Z_cor<-array(NA,c(n.run,2,4,10,2))
dim(TypeI)
m<-5000            #number of features 
n<-100            #sample size 
nx<-1           # number of x in simulation
n.hidden<-1   # number of z in simulation
lenx2<-c(0,0.5)*n  #two level of dependences between x1 and x2     
#meanz1<-c(0,0.4,0.8,1.2) #four levels of dependences between x1 and z1 
xz_cor<-c(0.1,0.3,0.6,0.9)
iter<-1;k<-1;h<-4
##track 
Bnn<-array(NA,c(n.run,2,4))
Bnn2<-array(NA,c(n.run,2,4))
Bnn3<-array(NA,c(n.run,2,4))
######################################
for (iter in 1:n.run){if(iter %%1==0) {print(10000+iter)}
  for (k in 1:1){
    for (h in 1:4){
      set.seed(1234+iter+k+h)
      ######################0.dataset####################
      e<-matrix(rnorm(n*m),n,m)
      dat=mvrnorm(n = n, rep(0, 2), Sigma=matrix(c(1,xz_cor[h],xz_cor[h],1),2,2))
      x=as.matrix(dat[,1])
      z=as.matrix(dat[,2])
      beta1<-matrix(0,m,1)
      ##delta matrix: coeffs matrix for z 
      #2000 features are influenced by the hidden variables 
      delta<-matrix(0,m,1)
      Sigma <- matrix(c(1,0,0,1),2,2)   #B D correlation
      coef<-mvrnorm(n=m,mu=c(0,0),Sigma=Sigma)
      cor(coef[,1],coef[,2])
      index1<-sample(c(1:m),size=2000,replace=F)  #sparsity
      index2<-sample(setdiff(c(1:m),index1),size=200,replace=F)
      index3<-sample(c(index1),size=1800,replace=F)
      beta1[index1]<-coef[index1,1]
      delta[index3,1]<-coef[index3,2]
      delta[index2,1]<-coef[index2,2]
      beta=beta1
      y<-matrix(NA,n,m)
      y<-x%*%(t(beta))+z%*%(t(delta))+e
      
      ########################1.Known#######################################################
      buja.num <- num.sv(t(y),model.matrix(~1+x), B =20, method = "be", seed = NULL)
      Bnn[iter,k,h]<-buja.num
      buja.over<-nx+n.hidden
      ########
      measure<-Get.measure(Y=y,X=x,Hidden=z,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,1,1]<-measure$mse.output)  
      (NUM[iter,k,h,1,1]<-measure$num.output )
      (NUMT[iter,k,h,1,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,1,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,1,1]<-measure$auc.output 
      ###power
      power[iter,k,h,1,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,1,1]<-measure$type1.output
      #####overestimated#####
      measure<-Get.measure(Y=y,X=x,Hidden=z,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,1,2]<-measure$mse.output)  
      (NUM[iter,k,h,1,2]<-measure$num.output )
      (NUMT[iter,k,h,1,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,1,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,1,2]<-measure$auc.output 
      ###power
      power[iter,k,h,1,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,1,2]<-measure$type1.output
      ########################2.NotAdj#######################################################
      measure<-Get.measure(Y=y,X=x,Hidden=NULL,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,2,1]<-measure$mse.output)  
      (NUM[iter,k,h,2,1]<-measure$num.output )
      (NUMT[iter,k,h,2,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,2,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,2,1]<-measure$auc.output 
      ###power
      power[iter,k,h,2,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,2,1]<-measure$type1.output
      
      #####overestimated#####
      measure<-Get.measure(Y=y,X=x,Hidden=z,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,2,2]<-measure$mse.output)  
      (NUM[iter,k,h,2,2]<-measure$num.output )
      (NUMT[iter,k,h,2,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,2,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,2,2]<-measure$auc.output 
      ###power
      power[iter,k,h,2,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,2,2]<-measure$type1.output
      ########################3.Two-setp SVA#######################################################
      sva.result<-sva(dat=t(y),mod=model.matrix(~x),method ="two-step",n.sv=buja.num)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=sva.result$sv,Coef=beta1) 
      Z_the[iter,k,h,1,1]=abs(sum(z*sva.result$sv))/(Norm(z,2)*Norm(sva.result$sv,2))
      Z_cor[iter,k,h,1,1]=summary(lm(z~sva.result$sv))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,3,1]<-measure$mse.output)  
      (NUM[iter,k,h,3,1]<-measure$num.output )
      (NUMT[iter,k,h,3,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,3,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,3,1]<-measure$auc.output 
      ###power
      power[iter,k,h,3,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,3,1]<-measure$type1.output
      
      #####overestimated#####
      sva.result<-sva(dat=t(y),mod=model.matrix(~x),method ="two-step",n.sv= buja.over)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=sva.result$sv,Coef=beta1) 
      #Z_the[iter,k,h,1,2]=abs(sum(z*sva.result$sv))
      Z_cor[iter,k,h,1,2]=summary(lm(z~sva.result$sv))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,3,2]<-measure$mse.output)  
      (NUM[iter,k,h,3,2]<-measure$num.output )
      (NUMT[iter,k,h,3,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,3,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,3,2]<-measure$auc.output 
      ###power
      power[iter,k,h,3,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,3,2]<-measure$type1.output
      ##########################4. Iter SVA #################################
      isva.result<-sva(dat=t(y),mod=model.matrix(~x),method ="irw",n.sv =buja.num)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=isva.result$sv,Coef=beta1) 
      Z_the[iter,k,h,2,1]=abs(sum(z*isva.result$sv))/(Norm(z,2)*Norm(isva.result$sv,2))
      Z_cor[iter,k,h,2,1]=summary(lm(z~isva.result$sv))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,4,1]<-measure$mse.output)  
      (NUM[iter,k,h,4,1]<-measure$num.output )
      (NUMT[iter,k,h,4,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,4,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,4,1]<-measure$auc.output 
      ###power
      power[iter,k,h,4,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,4,1]<-measure$type1.output
      
      #####overestimated#####
      isva.result<-sva(dat=t(y),mod=model.matrix(~x),method ="irw",n.sv =buja.over)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=isva.result$sv,Coef=beta1) 
      #Z_the[iter,k,h,2,2]=abs(sum(z*isva.result$sv))
      Z_cor[iter,k,h,2,2]=summary(lm(z~isva.result$sv))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,4,2]<-measure$mse.output)  
      (NUM[iter,k,h,4,2]<-measure$num.output )
      (NUMT[iter,k,h,4,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,4,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,4,2]<-measure$auc.output 
      ###power
      power[iter,k,h,4,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,4,2]<-measure$type1.output
      
      ####################5.PCA-buja #####################################
      n.pca<-buja.num
      ypca<-prcomp(y, retx = TRUE, center = T, scale. = T, tol = NULL)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=ypca$x[,1:as.numeric(n.pca)],Coef=beta1) 
      Z_the[iter,k,h,3,1]=abs(sum(z*ypca$x[,1:as.numeric(n.pca)]))/(Norm(z,2)*Norm(ypca$x[,1:as.numeric(n.pca)],2))
      Z_cor[iter,k,h,3,1]=summary(lm(z~ypca$x[,1:as.numeric(n.pca)]))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,5,1]<-measure$mse.output)  
      (NUM[iter,k,h,5,1]<-measure$num.output )
      (NUMT[iter,k,h,5,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,5,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,5,1]<-measure$auc.output 
      ###power
      power[iter,k,h,5,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,5,1]<-measure$type1.output
      
      #####overestimated#####
      n.pca<-buja.over
      measure<-Get.measure(Y=y,X=x,Hidden=ypca$x[,1:as.numeric(n.pca)],Coef=beta1) 
      #Z_the[iter,k,h,3,1]=abs(sum(z*ypca$x[,1:as.numeric(n.pca)]))
      Z_cor[iter,k,h,3,2]=summary(lm(z~ypca$x[,1:as.numeric(n.pca)]))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,5,2]<-measure$mse.output)  
      (NUM[iter,k,h,5,2]<-measure$num.output )
      (NUMT[iter,k,h,5,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,5,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,5,2]<-measure$auc.output 
      ###power
      power[iter,k,h,5,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,5,2]<-measure$type1.output
      
      ############################6. rPCA ################################
      y.residual<-y-lm.mp(y, formula=~x, store.fitted = T)$fitted
      y.residual.pca<-prcomp(y.residual, retx = TRUE, center = T, scale. = T, tol = NULL)
      ####
      Z_the[iter,k,h,4,1]=abs(sum(z*y.residual.pca$x[,1:as.numeric(buja.num)]))/(Norm(z,2)*Norm(y.residual.pca$x[,1:as.numeric(buja.num)],2))
      Z_cor[iter,k,h,4,1]=summary(lm(z~y.residual.pca$x[,1:as.numeric(buja.num)]))$r.squared
      measure<-Get.measure(Y=y,X=x,Hidden=y.residual.pca$x[,1:as.numeric(buja.num)],Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,6,1]<-measure$mse.output)  
      (NUM[iter,k,h,6,1]<-measure$num.output )
      (NUMT[iter,k,h,6,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,6,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,6,1]<-measure$auc.output 
      ###power
      power[iter,k,h,6,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,6,1]<-measure$type1.output
      
      #####overestimated#####
      measure<-Get.measure(Y=y,X=x,Hidden=y.residual.pca$x[,1:as.numeric(buja.over)],Coef=beta1) 
      ###MSE measure 
      #Z_the[iter,k,h,4,1]=abs(sum(z*y.residual.pca$x[,1:as.numeric(buja.num)]))
      Z_cor[iter,k,h,4,2]=summary(lm(z~y.residual.pca$x[,1:as.numeric(buja.over)]))$r.squared
      (MSE[iter,k,h,6,2]<-measure$mse.output)  
      (NUM[iter,k,h,6,2]<-measure$num.output )
      (NUMT[iter,k,h,6,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,6,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,6,2]<-measure$auc.output 
      ###power
      power[iter,k,h,6,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,6,2]<-measure$type1.output
      
      ##############################7. LFMM #############################
      
      #dSVA.result<-Get_SV_New(Y=y,X=x, ncomp=buja.num)
      lfmm.result<-lfmm(y, x, k=buja.num)
      ####
      Z_the[iter,k,h,5,1]=abs(sum(z*lfmm.result))/(Norm(z,2)*Norm(lfmm.result,2))
      Z_cor[iter,k,h,5,1]=summary(lm(z~lfmm.result))$r.squared
      measure<-Get.measure(Y=y,X=x,Hidden=lfmm.result,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,7,1]<-measure$mse.output)  
      (NUM[iter,k,h,7,1]<-measure$num.output )
      (NUMT[iter,k,h,7,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,7,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,7,1]<-measure$auc.output 
      ###power
      power[iter,k,h,7,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,7,1]<-measure$type1.output
      
      
      #####overestimated#####
      #dSVA.result<-Get_SV_New(Y=y,X=x, ncomp=buja.over)
      LFMM.result<-lfmm(y, x, k=buja.over)
      ####
      measure<-Get.measure(Y=y,X=x,Hidden=LFMM.result,Coef=beta1)
      #Z_the[iter,k,h,5,1]=abs(sum(z*lfmm.result))
      Z_cor[iter,k,h,5,2]=summary(lm(z~lfmm.result))$r.squared
      ###MSE measure 
      (MSE[iter,k,h,7,2]<-measure$mse.output)  
      (NUM[iter,k,h,7,2]<-measure$num.output )
      (NUMT[iter,k,h,7,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,7,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,7,2]<-measure$auc.output 
      ###power
      power[iter,k,h,7,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,7,2]<-measure$type1.output
      
      
      #############################8. RUV4-ncol #############################
      length((which(beta1!=0&delta[,1]!=0)))
      ctl0<-c(sample(which(beta1!=0&delta[,1]!=0),size=100,replace=F),
              sample(which(beta1==0&delta[,1]!=0),size=100,replace=F),
              sample(which(beta1==0&delta[,1]==0),size=100,replace=F),
              sample(which(beta1!=0&delta[,1]==0),size=100,replace=F))
      ctl<-rep(F,m)
      ctl[ctl0]<-T;sum(ctl)
      #
      #(opt.k<-getK(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=as.matrix(x[,2]))$k)
      (opt.k<-getK(Y=y,X=as.matrix(x),ctl=ctl)$k)
      #opt.k<-2
      #ruv4.result<-RUV4(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=as.matrix(x[,2]),k=opt.k)
      ruv4.result<-RUV4(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=1,k=opt.k)
      if (opt.k>0)
      {
        Z_the[iter,k,h,6,1]=abs(sum(z*ruv4.result$W))/(Norm(z,2)*Norm(ruv4.result$W,2))
        Z_cor[iter,k,h,6,1]=summary(lm(z~ruv4.result$W))$r.squared
      }
      ####
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(ruv4.result$betahat),Coef=as.vector(beta1),Pvalue=ruv4.result$p) 
      ###MSE measure 
      (MSE[iter,k,h,8,1]<-measure$mse.output)  
      (NUM[iter,k,h,8,1]<-measure$num.output )
      (NUMT[iter,k,h,8,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,8,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,8,1]<-measure$auc.output 
      ###power
      power[iter,k,h,8,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,8,1]<-measure$type1.output
      
      #####overestimated#####
      ruv4.result<-RUV4(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=1,k=(opt.k*2))
      ####
      #Z_the[iter,k,h,6,1]=abs(sum(z*ruv4.result$W))
      if (opt.k>0)
      {
        Z_cor[iter,k,h,6,2]=summary(lm(z~ruv4.result$W))$r.squared
      }
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(ruv4.result$betahat),Coef=as.vector(beta1),Pvalue=ruv4.result$p) 
      ###MSE measure 
      (MSE[iter,k,h,8,2]<-measure$mse.output)  
      (NUM[iter,k,h,8,2]<-measure$num.output )
      (NUMT[iter,k,h,8,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,8,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,8,2]<-measure$auc.output 
      ###power
      power[iter,k,h,8,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,8,2]<-measure$type1.output
      #############################9. RUV4-buja #############################
      opt.k<-buja.num
      ruv4.result<-RUV4(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=1,k=opt.k)
      ####
      Z_the[iter,k,h,7,1]=abs(sum(z*ruv4.result$W))/(Norm(z,2)*Norm(ruv4.result$W,2))
      Z_cor[iter,k,h,7,1]=summary(lm(z~ruv4.result$W))$r.squared
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(ruv4.result$betahat),Coef=as.vector(beta1),Pvalue=ruv4.result$p) 
      ###MSE measure 
      (MSE[iter,k,h,9,1]<-measure$mse.output)  
      (NUM[iter,k,h,9,1]<-measure$num.output )
      (NUMT[iter,k,h,9,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,9,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,9,1]<-measure$auc.output 
      ###power
      power[iter,k,h,9,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,9,1]<-measure$type1.output
      
      #####overestimated#####
      opt.k<-buja.over
      ruv4.result<-RUV4(Y=y,X=as.matrix(x[,1]),ctl=ctl,Z=1,k=(opt.k))
      ####
      #Z_the[iter,k,h,7,2]=abs(sum(z*ruv4.result$W))
      Z_cor[iter,k,h,7,2]=summary(lm(z~ruv4.result$W))$r.squared
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(ruv4.result$betahat),Coef=as.vector(beta1),Pvalue=ruv4.result$p) 
      ###MSE measure 
      (MSE[iter,k,h,9,2]<-measure$mse.output)  
      (NUM[iter,k,h,9,2]<-measure$num.output )
      (NUMT[iter,k,h,9,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,9,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,9,2]<-measure$auc.output 
      ###power
      power[iter,k,h,9,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,9,2]<-measure$type1.output
      #############################10. Leapp-RR #############################
      X.data <- data.frame(x = x)
      rr.result<-cate(formula=~x|1, X.data,Y=y, r=buja.num, 
                      fa.method ="ml", adj.method = "rr", 
                      psi = psi.huber, calibrate = F)
      ##
      Z_the[iter,k,h,8,1]=abs(sum(z*rr.result$Z))/(Norm(z,2)*Norm(rr.result$Z,2))
      Z_cor[iter,k,h,8,1]=summary(lm(z~rr.result$Z))$r.squared
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(rr.result$beta),Coef=as.vector(beta1),Pvalue=rr.result$beta.p.value) 
      ###MSE measure 
      (MSE[iter,k,h,10,1]<-measure$mse.output)  
      (NUM[iter,k,h,10,1]<-measure$num.output )
      (NUMT[iter,k,h,10,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,10,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,10,1]<-measure$auc.output 
      ###power
      power[iter,k,h,10,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,10,1]<-measure$type1.output
      
      
      #####overestimated#####
      X.data <- data.frame(x = x)
      rr.result<-cate(formula=~x|1, X.data, Y=y, r=buja.over, 
                      fa.method ="ml", adj.method = "rr", 
                      psi = psi.huber, calibrate = F)
      ##
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(rr.result$beta),Coef=as.vector(beta1),Pvalue=rr.result$beta.p.value) 
      ###MSE measure 
      #Z_the[iter,k,h,8,1]=abs(sum(z*rr.result$Z))
      Z_cor[iter,k,h,8,2]=summary(lm(z~rr.result$Z))$r.squared
      (MSE[iter,k,h,10,2]<-measure$mse.output)  
      (NUM[iter,k,h,10,2]<-measure$num.output )
      (NUMT[iter,k,h,10,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,10,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,10,2]<-measure$auc.output 
      ###power
      power[iter,k,h,10,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,10,2]<-measure$type1.output
      
      #############################11. Leapp-RR-DAM #############################
      mad.result<-cate(formula=~x|1, X.data,Y=y, r=buja.num, 
                       fa.method ="ml", adj.method = "rr", 
                       psi = psi.huber, calibrate = TRUE)
      
      ###
      Z_the[iter,k,h,9,1]=abs(sum(z*mad.result$Z))/(Norm(z,2)*Norm(mad.result$Z,2))
      Z_cor[iter,k,h,9,1]=summary(lm(z~mad.result$Z))$r.squared
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(mad.result$beta),Coef=as.vector(beta1),Pvalue=mad.result$beta.p.value) 
      ###MSE measure 
      (MSE[iter,k,h,11,1]<-measure$mse.output)  
      (NUM[iter,k,h,11,1]<-measure$num.output )
      (NUMT[iter,k,h,11,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,11,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,11,1]<-measure$auc.output 
      ###power
      power[iter,k,h,11,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,11,1]<-measure$type1.output
      
      #####overestimated#####
      mad.result<-cate(formula=~x|1, X.data,Y=y, r=buja.over, 
                       fa.method ="ml", adj.method = "rr", 
                       psi = psi.huber, calibrate = TRUE)
      ###
      #Z_the[iter,k,h,9,1]=abs(sum(z*mad.result$Z))
      Z_cor[iter,k,h,9,2]=summary(lm(z~mad.result$Z))$r.squared
      measure<-Get.measure2(Y=y,X=x,Estimator=as.vector(mad.result$beta),Coef=as.vector(beta1),Pvalue=mad.result$beta.p.value) 
      ###MSE measure 
      (MSE[iter,k,h,11,2]<-measure$mse.output)  
      (NUM[iter,k,h,11,2]<-measure$num.output )
      (NUMT[iter,k,h,11,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,11,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,11,2]<-measure$auc.output 
      ###power
      power[iter,k,h,11,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,11,2]<-measure$type1.output
      ##############################12. dSVA #############################
      #dSVA.result<-Get_SV_New(Y=y,X=x, ncomp=buja.num)
      dSVA.result<-dSVA(y, x, ncomp=buja.num)
      ####
      Z_the[iter,k,h,10,1]=abs(sum(z*dSVA.result$Z))/(Norm(z,2)*Norm(dSVA.result$Z,2))
      Z_cor[iter,k,h,10,1]=summary(lm(z~dSVA.result$Z))$r.squared
      measure<-Get.measure(Y=y,X=x,Hidden=dSVA.result$Z,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,12,1]<-measure$mse.output)  
      (NUM[iter,k,h,12,1]<-measure$num.output )
      (NUMT[iter,k,h,12,1]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,12,1]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,12,1]<-measure$auc.output 
      ###power
      power[iter,k,h,12,1]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,12,1]<-measure$type1.output
      
      
      #####overestimated#####
      #dSVA.result<-Get_SV_New(Y=y,X=x, ncomp=buja.over)
      dSVA.result<-dSVA(y, x, ncomp=buja.over)
      ####
      #Z_the[iter,k,h,10,1]=abs(sum(z*dSVA.result$Z))
      Z_cor[iter,k,h,10,2]=summary(lm(z~dSVA.result$Z))$r.squared
      measure<-Get.measure(Y=y,X=x,Hidden=dSVA.result$Z,Coef=beta1) 
      ###MSE measure 
      (MSE[iter,k,h,12,2]<-measure$mse.output)  
      (NUM[iter,k,h,12,2]<-measure$num.output )
      (NUMT[iter,k,h,12,2]<-measure$num.true.output)
      ###FDR measure 
      (FDR[iter,k,h,12,2]<-measure$fdr.output)
      ###AUC measure 
      AUC[iter,k,h,12,2]<-measure$auc.output 
      ###power
      power[iter,k,h,12,2]<-measure$power.output
      ###Type I error
      TypeI[iter,k,h,12,2]<-measure$type1.output
      
      
      
    }
  }
}






save.image("Data-Case9.RData")

