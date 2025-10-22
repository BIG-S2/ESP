library(Matrix)
library(vows)
library(sva)
library(AUC)
library(ruv) 
##
library(pscl)#for inverse gamma 
library(MASS)#for normal  
##
#library(jackstraw)
#library(qvalue)
library(xtable)
##
library(cate)
##
library(bindata)
##
library(lfmm)
#library(iterSVA)
library(dSVA)
########################################################################
########################performances####################################
########################################################################
###MSE measure, FDR measure, AUC measure, Power measure  
Get.measure<-function(Y,X,Hidden,Coef){
  if(is.null(Hidden)){
    Coef.hat<-lm.mp(Y, formula=~X, store.fitted = FALSE)$coef[2,]
    pvalue.raw<-F.mp(Y~X,which=2)$pvalue
  }
  if(!is.null(Hidden)){
    Coef.hat<-lm.mp(Y, formula=~X+Hidden, store.fitted = FALSE)$coef[2,]
    pvalue.raw<-F.mp(Y~X+Hidden,which=2)$pvalue
  }
  ##MSE
  mse.output<-sum((Coef.hat-Coef)^2)/length(Coef)  
  ##FDR
  pvalue.adjust<-p.adjust(pvalue.raw, method = "fdr")
  fdr.output<-sum((Coef==0&pvalue.adjust<=0.05)*1)/sum((pvalue.adjust<=0.05)*1) 
  if (sum((pvalue.adjust<0.05)*1)==0){ 
    fdr.output<-0
  }
  num.output<-sum((pvalue.adjust<=0.05)*1) 
  num.true.output<-sum((Coef!=0&pvalue.adjust<=0.05)*1)
  ###
  # mse_sub1.output<-NA
  # mse_sub2.output<-NA
  # if(sum((pvalue.adjust<0.05)*1)>0){
  #   false.pos<-(Coef==0)&(pvalue.adjust<=0.05)
  #   ##MSE-sub1
  #   mse_sub1.output<-sum((Coef.hat[false.pos]-Coef[false.pos])^2)/length(Coef[false.pos])  
  #   ##MSE-sub2
  #   mse_sub2.output<-sum((Coef.hat[!false.pos]-Coef[!false.pos])^2)/length(Coef[!false.pos])  
  #   
  # }
  ##AUC
  labels<-(Coef!=0)*1
  predictions<-(pvalue.raw<0.05)*1
  raw.output<-sum(predictions)
  auc.output<-AUC::auc(roc(predictions,as.factor(labels))) 
  #Type I error
  type1.output<-sum((pvalue.raw<0.05&Coef==0)*1)/length(Coef) 
  ##Power
  power.output<-sum((pvalue.raw<0.05&Coef!=0)*1)/sum(Coef!=0)
  ###
  return(list(mse.output=mse.output,
              fdr.output=fdr.output,
              num.output=num.output,
              raw.output=raw.output,
              auc.output=auc.output,
              type1.output=type1.output,
              num.true.output=num.true.output,
              power.output=power.output,
              Coef.hat=Coef.hat))
}

####used for other methods: RUV,LEAPP-RR
Get.measure2<-function(Y,X,Estimator,Coef,Pvalue){
  
  ##MSE
  mse.output<-sum((Estimator-Coef)^2)/length(Coef)  
  ##FDR
  pvalue.raw<-Pvalue
  pvalue.adjust<-p.adjust(pvalue.raw, method = "fdr")
  fdr.output<-sum((Coef==0&pvalue.adjust<=0.05)*1)/sum((pvalue.adjust<=0.05)*1) 
  if (sum((pvalue.adjust<0.05)*1)==0){ 
    fdr.output<-0}
  num.output<-sum((pvalue.adjust<=0.05)*1) 
  num.true.output<-sum((Coef!=0&pvalue.adjust<=0.05)*1)
  ##
  ##AUC
  labels<-(Coef!=0)*1
  predictions<-(pvalue.raw<0.05)*1
  raw.output<-sum(predictions)
  auc.output<-auc(roc(predictions,as.factor(labels))) 
  #Type I 
  type1.output<-sum((pvalue.raw<0.05&Coef==0)*1)/length(Coef) 
  ##
  ##Power
  power.output<-sum((pvalue.raw<0.05&Coef!=0)*1)/sum(Coef!=0)
  return(list(mse.output=mse.output,
              fdr.output=fdr.output,
              num.output=num.output,
              raw.output=raw.output,
              auc.output=auc.output,
              type1.output=type1.output,
              num.true.output=num.true.output,
              power.output=power.output))
}#######

###########################################################################
############################################################################
################################# Functions ################################
############################################################################
############################################################################


rPCA20_hidden<-function(Y,X,nComp,nIter,mCoverge,pLevelx,pLevelz){
  ###read in 
  Y<-as.matrix(Y);
  X<-as.matrix(X);
  n<-dim(Y)[1];m<-dim(Y)[2] 
  nx<-dim(X)[2]   ##number of covariates 
  n.rpca2<-nComp  ##number of top PCs taken as Z 
  n.iter<-nIter   ##repeated times of iterative step 
  pcutx<-pLevelx      ##significance level for x
  pcutz<-pLevelz      ##significance level for z
  pcutw<-0.05         ##significance level for weight vector 
  convcut<-mCoverge  ##converge cut for Y2
  if(n.rpca2==0){
    out<-NULL
  } 
  if(n.rpca2>=1){
    ################ Get Z1.est by All Features ################
    ####A.Detect significant PCs for observed covariates  
    ypca<-prcomp(Y, retx = TRUE, center = T, scale. = T, tol = NULL)
    pvalue.raw.x<-matrix(NA,(min(n,m)*0.5),nx)
    for ( ii in 1:(min(n,m)*0.5)){
      for( jj in 1:nx){
        pvalue.raw.x[ii,jj]<-summary(lm(X[,jj]~ypca$x[,ii]))$coef[2,4]
      }
    }
    pvalue.adj.x<-apply(pvalue.raw.x,2,function(xx){p.adjust(xx, method = "fdr")})
    sign.pc<-unique(which(pvalue.adj.x<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
    print(paste("On Y, selected PCs for X:"))
    print(sign.pc)
    num.pc<-length(sign.pc)
    if(num.pc==0){
      print('On Y, number of selected PC is 0, will take top PC as Z')
      Z1.est<-ypca$x[,1:n.rpca2]
    }
    if(num.pc!=0){
      ####B.Find the weight vector, suggested by Seunggeun
      ##Get p-value in model 1 and model 2 
      ##Get the weight based on the two p-values 
      pvalue.y.m1<-matrix(NA,num.pc,m)
      pvalue.y.m2<-matrix(NA,num.pc,m)
      weight<-matrix(NA,num.pc,m)
      for (kk in 1: num.pc){
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]],which=2)$pvalue
        pvalue.y.m1[kk,]<-p.adjust(temp, method = "fdr")
        ##
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]]+X,which=2)$pvalue
        pvalue.y.m2[kk,]<-p.adjust(temp, method = "fdr")
        ##
        weight[kk,]<-0
        #weight[kk,]<-(ypca$rot[,sign.pc[kk]])^2
        weight[kk,(pvalue.y.m1[kk,]<=pcutw)&(pvalue.y.m2[kk,]>=pcutw)]<-1
      }
      ####C. rPCA2 based on adjusted Beta
      ## Adjust beta by weight vector 
      beta.ori<-lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[c(2:(num.pc+1)),]
      beta.adj<-rbind(lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[1,],beta.ori*weight)
      ## Get the residual PCs
      yres<-Y-cbind(1,ypca$x[,sign.pc])%*%beta.adj
      yres.pca<-prcomp(yres, retx = TRUE, center = T, scale = T, tol = NULL)
      Z1.est<-yres.pca$x[,1:n.rpca2]
    }
    ################ Update Z1.est by Y2 (Subset of Features) ##############
    TempZ<-Z1.est
    for (gg in 1:n.iter){
      ####D. Get Y2: subset of Y with large Z and low X 
      pvalue.y2<-matrix(NA,2,m)
      temp<-F.mp(Y~X+TempZ,which=c(2:(1+nx)))$pvalue
      pvalue.y2[1,]<-p.adjust(temp, method = "fdr")
      temp<-F.mp(Y~TempZ+X,which=c(2:(1+n.rpca2)))$pvalue
      pvalue.y2[2,]<-p.adjust(temp, method = "fdr")
      temp<-(pvalue.y2[1,]>= pcutx)&(pvalue.y2[2,]<=pcutz)
      if(gg>1){
        diff.y2<-sum((temp.y2*1)!=(temp*1))
        print(paste("Number of different features between new Y2 and old Y2: ",diff.y2));
        if(diff.y2<=convcut){
          print('Y2 Converged, Done')
          break;
        }
      }
      temp.y2<-temp
      Y2<-Y[,temp]
      m2<-dim(Y2)[2]
      print(paste("Iter on Y2:",gg,", Dim :", m2,", n.hiden :", n.rpca2))
      ####E. Detect significant PCs for observed covariates on Y2
      ypca2<-prcomp(Y2, retx = TRUE, center = T, scale. = T, tol = NULL)
      pvalue.raw.x2<-matrix(NA,(n*0.5),nx)
      for ( ii in 1:(min(n,m2)*0.5)){
        for( jj in 1:nx){
          pvalue.raw.x2[ii,jj]<-summary(lm(X[,jj]~ypca2$x[,ii]))$coef[2,4]
        }
      }
      pvalue.adj.x2<-apply(pvalue.raw.x2,2,function(xx){p.adjust(xx, method = "fdr")})
      sign.pc2<-unique(which(pvalue.adj.x2<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
      print('The PCs related to X in Y2: ');
      print(sign.pc2)
      num.pc2<-length(sign.pc2)
      ## If no selected PC, output the PCs
      if(num.pc2==0){
        print('On Y2, number of selected PC is 0, will take top PC as Z')
        Z2.est<-ypca2$x[,1:n.rpca2]
        TempZ<-Z2.est
      }
      ###If there is selected PC, do rPCA2 again  
      if(num.pc2>0){
        ## Get p-value in model 1 and model 2 
        ## Get the weight based on the two p-values 
        pvalue.y2.m1<-matrix(NA,num.pc2,m2)
        pvalue.y2.m2<-matrix(NA,num.pc2,m2)
        weight2<-matrix(NA,num.pc2,m2)
        for (kk in 1: num.pc2){
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]],which=2)$pvalue
          pvalue.y2.m1[kk,]<-p.adjust(temp, method = "fdr")
          ##
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]]+X,which=2)$pvalue
          pvalue.y2.m2[kk,]<-p.adjust(temp, method = "fdr")
          ##
          weight2[kk,]<-0
          #weight2[kk,]<-(ypca2$rot[,sign.pc2[kk]])^2
          weight2[kk,(pvalue.y2.m1[kk,]<=pcutw)&(pvalue.y2.m2[kk,]>=pcutw)]<-1
        }
        print(paste('Number of features with at least one weight2==1:',sum(apply(weight2,2,sum)>0)));
        ####F. rPCA2 based on adjusted Beta
        ##Adjust the Beta by weight vector 
        beta.ori2<-lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[c(2:(num.pc2+1)),]
        beta.adj2<-rbind(lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[1,],beta.ori2*weight2)
        ##Get the residual PCs
        yres2<-Y2-cbind(1,ypca2$x[,sign.pc2])%*%beta.adj2
        yres.pca2<-prcomp(yres2, retx = TRUE, center = T, scale = T, tol = NULL)
        Z2.est<-yres.pca2$x[,1:c(n.rpca2)]
        TempZ<-Z2.est
      }
    }
    out<-list(Z1.est=Z1.est,Z2.est=Z2.est,temp.y2=temp.y2)
  }
  return(out)
}



#########################################################################
########################################################################
#########################################################################
########################################################################


rPCA29_hidden<-function(Y,X,nComp,nIter,mCoverge,pLevelx,pLevelz,pselect=0.2,nDetect=15){
  ###read in 
  Y<-as.matrix(Y);
  X<-as.matrix(X);
  n<-dim(Y)[1];m<-dim(Y)[2] 
  #xc<-Covariates #observed covariates
  nx<-dim(X)[2]   ##number of covariates 
  n.rpca2<-nComp  ##number of top PCs taken as Z 
  n.iter<-nIter   ##repeated times of iterative step 
  pcutx<-pLevelx      ##significance level for x
  pcutz<-pLevelz      ##significance level for z
  pcutw<-0.05        ##significance level for weight vector 
  pcuthidden<-pselect    ##detect remaining surrgoate of X 
  convcut<-mCoverge  ##converge cut for Y2
  if(n.rpca2==0){
    out<-NULL
  } 
  if(n.rpca2>=1){
    ################ Get Z1.est by All Features ################
    ####A.Detect significant PCs for observed covariates  
    ypca<-prcomp(Y, retx = TRUE, center = T, scale. = T, tol = NULL)
    pvalue.raw.x<-matrix(NA,(min(n,m)*0.5),nx)
    for ( ii in 1:(min(n,m)*0.5)){
      for( jj in 1:nx){
        pvalue.raw.x[ii,jj]<-summary(lm(X[,jj]~ypca$x[,ii]))$coef[2,4]
      }
    }
    pvalue.adj.x<-apply(pvalue.raw.x,2,function(xx){p.adjust(xx, method = "fdr")})
    sign.pc<-unique(which(pvalue.adj.x<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
    print(paste("On Y, selected PCs for X:"))
    print(sign.pc)
    num.pc<-length(sign.pc)
    if(num.pc==0){
      print('On Y, number of selected PC is 0, will take top PCs as Z')
      Z1.est<-ypca$x[,1:n.rpca2]   
    }
    if(num.pc!=0){
      ####B.Find the weight vector, suggested by Seunggeun
      ##Get p-value in model 1 and model 2; Get the weight based on the two p-values 
      pvalue.y.m1<-matrix(NA,num.pc,m)
      pvalue.y.m2<-matrix(NA,num.pc,m)
      weight<-matrix(NA,num.pc,m)
      for (kk in 1: num.pc){
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]],which=2)$pvalue
        pvalue.y.m1[kk,]<-p.adjust(temp, method = "fdr")
        ##
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]]+X,which=2)$pvalue
        pvalue.y.m2[kk,]<-p.adjust(temp, method = "fdr")
        ##
        weight[kk,]<-0
        weight[kk,(pvalue.y.m1[kk,]<=pcutw)&(pvalue.y.m2[kk,]>=pcutw)]<-1
      }
      ####C. rPCA2 based on adjusted Beta
      ## Adjust beta by weight vector 
      beta.ori<-lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[c(2:(num.pc+1)),]
      beta.adj<-rbind(lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[1,],beta.ori*weight)
      ## Get the residual PCs
      yres<-Y-cbind(1,ypca$x[,sign.pc])%*%beta.adj
      yres.pca<-prcomp(yres, retx = TRUE, center = T, scale = T, tol = NULL)
      Z1.est<-yres.pca$x[,1:n.rpca2]
      ####
    }
    ################ Update Z1.est by Y2 (Subset of Features) ##############
    TempZ<-Z1.est
    ###Do iteration 
    for (gg in 1:n.iter){
      ####D. Get Y2: subset of Y with ZD>>XB
      pvalue.y2<-matrix(NA,2,m)
      temp<-F.mp(Y~X+TempZ,which=c(2:(1+nx)))$pvalue
      pvalue.y2[1,]<-p.adjust(temp, method = "fdr")
      temp<-F.mp(Y~TempZ+X,which=c(2:(1+n.rpca2)))$pvalue
      pvalue.y2[2,]<-p.adjust(temp, method = "fdr")
      temp<-(pvalue.y2[1,]>= pcutx)&(pvalue.y2[2,]<=pcutz)
      if(gg>1){
        diff.y2<-sum((temp.y2*1)!=(temp*1))
        print(paste("Difference (new Y2 and old Y2): ",diff.y2));
        if(diff.y2<=convcut){
          print('Y2 Converged, Done')
          break;
        }
      }
      temp.y2<-temp
      Y2<-Y[,temp]
      m2<-dim(Y2)[2]
      print(paste("Iter on Y2:",gg,", Dim :", m2,", n.hidden :", n.rpca2))
      ####E. Detect significant PCs for observed covariates on Y2
      ypca2<-prcomp(Y2, retx = TRUE, center = T, scale. = T, tol = NULL)
      pvalue.raw.x2<-matrix(NA,(min(n,m2)*0.5),nx)
      for ( ii in 1:(min(n,m2)*0.5)){
        for( jj in 1:nx){
          pvalue.raw.x2[ii,jj]<-summary(lm(X[,jj]~ypca2$x[,ii]))$coef[2,4]
        }
      }
      pvalue.adj.x2<-apply(pvalue.raw.x2,2,function(xx){p.adjust(xx, method = "fdr")})
      sign.pc2<-unique(which(pvalue.adj.x2<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
      print(paste('PCs from Y2 related to X: '));print(sign.pc2)
      num.pc2<-length(sign.pc2)
      ## If no selected PC, output the PCs
      if(num.pc2==0){
        print('On Y2, number of selected PC is 0, will take top PCs as Z')
        yres2<-Y2
        yres.pca2<-ypca2
        Z2.est0<-ypca2$x[,1:n.rpca2]
        TempZ<-Z2.est0
      }
      ###If there is selected PC, do rPCA2 again  
      if(num.pc2>0){
        ## Get p-value in model 1 and model 2 
        ## Get the weight based on the two p-values 
        pvalue.y2.m1<-matrix(NA,num.pc2,m2)
        pvalue.y2.m2<-matrix(NA,num.pc2,m2)
        weight2<-matrix(NA,num.pc2,m2)
        for (kk in 1: num.pc2){
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]],which=2)$pvalue
          pvalue.y2.m1[kk,]<-p.adjust(temp, method = "fdr")
          ##
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]]+X,which=2)$pvalue
          pvalue.y2.m2[kk,]<-p.adjust(temp, method = "fdr")
          ##
          weight2[kk,]<-0
          weight2[kk,(pvalue.y2.m1[kk,]<=pcutw)&(pvalue.y2.m2[kk,]>=pcutw)]<-1
        }
        print(paste('Number of features with at least weight 1:',sum(apply(weight2,2,sum)>0)));
        ####F. rPCA2 based on adjusted Beta
        ##Adjust the Beta by weight vector 
        beta.ori2<-lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[c(2:(num.pc2+1)),]
        beta.adj2<-rbind(lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[1,],beta.ori2*weight2)
        ##Get the residual PCs
        yres2<-Y2-cbind(1,ypca2$x[,sign.pc2])%*%beta.adj2
        yres.pca2<-prcomp(yres2, retx = TRUE, center = T, scale = T, tol = NULL)
        Z2.est0<-yres.pca2$x[,1:c(n.rpca2)]
        ####
        TempZ<-Z2.est0
        ####
      }
    }
    ########Final Update on Weighted Residuals of Y2 ########
    n.test<-min(nDetect,min(n,m))
    pvalue.yres.n1<-matrix(NA,n.test,m2)
    pvalue.yres.n2<-matrix(NA,n.test,m2)
    part.hidden<-rep(NA,n.test)
    for (kk in 1: n.test){
      temp<-F.mp(Y2~yres.pca2$x[,kk],which=c(2))$pvalue
      pvalue.yres.n1[kk,]<-temp
      ##
      temp<-F.mp(Y2~yres.pca2$x[,kk]+X,which=c(2))$pvalue
      pvalue.yres.n2[kk,]<-temp
      ##
      if(sum((pvalue.yres.n1[kk,]<=pcutw))>10){
        part.hidden[kk]<-(sum((pvalue.yres.n1[kk,]<=pcutw)&(pvalue.yres.n2[kk,]>=pcutw))/sum((pvalue.yres.n1[kk,]<=pcutw)))
      }
      if(sum((pvalue.yres.n1[kk,]<=pcutw))<=10){#avoid wired case 
        part.hidden[kk]<-0
      }
    }
    print(part.hidden)
    diff.pc.r<-which(part.hidden>pcuthidden)
    if(length(diff.pc.r)>0){
      nt<-min(nx,length(diff.pc.r))
      remove.r<-sort(diff.pc.r,decreasing =T)[1:nt]
      ####
      print("Final update on Y2:");print(remove.r)
      yres3<-yres2-lm.mp(yres2, formula=~yres.pca2$x[, remove.r], store.fitted =T)$fitted
      yres.pca3<-prcomp(yres3, retx = TRUE, center = T, scale = T, tol = NULL)
      Z2.est<-yres.pca3$x[,1:c(n.rpca2)]
    } 
    if(length(diff.pc.r)==0){
      Z2.est<-Z2.est0
    }
    ########
    out<-list(Z1.est=Z1.est,Z2.est=Z2.est,temp.y2=temp.y2,detect.primary=part.hidden)
  }
  return(out)
}


#########################################################################
########################################################################
#########################################################################
########################################################################


rPCA30_hidden<-function(Y,X,nComp,nIter,mCoverge,pLevelx,pLevelz){
  ###read in 
  Y<-as.matrix(Y);
  X<-as.matrix(X);
  n<-dim(Y)[1];m<-dim(Y)[2] 
  nx<-dim(X)[2]   ##number of covariates 
  n.rpca2<-nComp  ##number of top PCs taken as Z 
  n.iter<-nIter   ##repeated times of iterative step 
  pcutx<-pLevelx      ##significance level for x
  pcutz<-pLevelz      ##significance level for z
  pcutw<-0.05        ##significance level for weight vector 
  convcut<-mCoverge  ##converge cut for Y2
  if(n.rpca2==0){
    out<-NULL
  } 
  if(n.rpca2>=1){
    ################ Get Z1.est by All Features ################
    ####A.Detect significant PCs for observed covariates  
    ypca<-prcomp(Y, retx = TRUE, center = T, scale. = T, tol = NULL)
    pvalue.raw.x<-matrix(NA,(min(n,m)*0.5),nx)
    for ( ii in 1:(min(n,m)*0.5)){
      for( jj in 1:nx){
        pvalue.raw.x[ii,jj]<-summary(lm(X[,jj]~ypca$x[,ii]))$coef[2,4]
      }
    }
    pvalue.adj.x<-apply(pvalue.raw.x,2,function(xx){p.adjust(xx, method = "fdr")})
    sign.pc<-unique(which(pvalue.adj.x<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
    print(paste("On Y, selected PCs for X:"))
    print(sign.pc)
    num.pc<-length(sign.pc)
    if(num.pc==0){
      print('On Y, number of selected PC is 0, will take top PCs as Z')
      Z1.est<-ypca$x[,1:n.rpca2]   
      yres<-Y
      yres.pca<-ypca
    }
    if(num.pc!=0){
      ####B.Find the weight vector, suggested by Seunggeun
      ##Get p-value in model 1 and model 2; Get the weight based on the two p-values 
      pvalue.y.m1<-matrix(NA,num.pc,m)
      pvalue.y.m2<-matrix(NA,num.pc,m)
      weight<-matrix(NA,num.pc,m)
      for (kk in 1: num.pc){
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]],which=2)$pvalue
        pvalue.y.m1[kk,]<-p.adjust(temp, method = "fdr")
        ##
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]]+X,which=2)$pvalue
        pvalue.y.m2[kk,]<-p.adjust(temp, method = "fdr")
        ##
        weight[kk,]<-0
        weight[kk,(pvalue.y.m1[kk,]<=pcutw)&(pvalue.y.m2[kk,]>=pcutw)]<-1
      }
      ####C. rPCA2 based on adjusted Beta
      ## Adjust beta by weight vector 
      beta.ori<-lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[c(2:(num.pc+1)),]
      beta.adj<-rbind(lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[1,],beta.ori*weight)
      ## Get the residual PCs
      yres<-Y-cbind(1,ypca$x[,sign.pc])%*%beta.adj
      yres.pca<-prcomp(yres, retx = TRUE, center = T, scale = T, tol = NULL)
      Z1.est<-yres.pca$x[,1:n.rpca2]
    }
    ################ Update Z1.est by Y2 (Subset of Features) ##############
    TempZ<-Z1.est
    ###Do iteration 
    for (gg in 1:n.iter){
      ####D. Get Y2: subset of Y with ZD>>XB
      pvalue.y2<-matrix(NA,2,m)
      temp<-F.mp(Y~X+TempZ,which=c(2:(1+nx)))$pvalue
      pvalue.y2[1,]<-p.adjust(temp, method = "fdr")
      temp<-F.mp(Y~TempZ+X,which=c(2:(1+n.rpca2)))$pvalue
      pvalue.y2[2,]<-p.adjust(temp, method = "fdr")
      temp<-(pvalue.y2[1,]>= pcutx)&(pvalue.y2[2,]<=pcutz)
      if(gg>1){
        diff.y2<-sum((temp.y2*1)!=(temp*1))
        print(paste("Difference (new Y2 and old Y2): ",diff.y2));
        if(diff.y2<=convcut){
          print('Y2 Converged, Done')
          break;
        }
      }
      temp.y2<-temp
      Y2<-Y[,temp]
      m2<-dim(Y2)[2]
      print(paste("Iter on Y2:",gg,", Dim :", m2,", n.hidden :", n.rpca2))
      ####E. Detect significant PCs for observed covariates on Y2
      ypca2<-prcomp(Y2, retx = TRUE, center = T, scale. = T, tol = NULL)
      pvalue.raw.x2<-matrix(NA,(min(n,m2)*0.5),nx)
      for ( ii in 1:(min(n,m2)*0.5)){
        for( jj in 1:nx){
          pvalue.raw.x2[ii,jj]<-summary(lm(X[,jj]~ypca2$x[,ii]))$coef[2,4]
        }
      }
      pvalue.adj.x2<-apply(pvalue.raw.x2,2,function(xx){p.adjust(xx, method = "fdr")})
      sign.pc2<-unique(which(pvalue.adj.x2<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
      print(paste('PCs from Y2 related to X: '));print(sign.pc2)
      num.pc2<-length(sign.pc2)
      ## If no selected PC, output the PCs
      if(num.pc2==0){
        print('On Y2, number of selected PC is 0, will take top PCs as Z')
        Z2.est<-ypca2$x[,1:n.rpca2]
        yres2<-Y2
        yres.pca2<-ypca2
      }
      ###If there is selected PC, do rPCA2 again  
      if(num.pc2>0){
        ## Get p-value in model 1 and model 2 
        ## Get the weight based on the two p-values 
        pvalue.y2.m1<-matrix(NA,num.pc2,m2)
        pvalue.y2.m2<-matrix(NA,num.pc2,m2)
        weight2<-matrix(NA,num.pc2,m2)
        for (kk in 1: num.pc2){
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]],which=2)$pvalue
          pvalue.y2.m1[kk,]<-p.adjust(temp, method = "fdr")
          ##
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]]+X,which=2)$pvalue
          pvalue.y2.m2[kk,]<-p.adjust(temp, method = "fdr")
          ##
          weight2[kk,]<-0
          weight2[kk,(pvalue.y2.m1[kk,]<=pcutw)&(pvalue.y2.m2[kk,]>=pcutw)]<-1
        }
        print(paste('Number of features removing any PC:',sum(apply(weight2,2,sum)>0)));
        ####F. rPCA2 based on adjusted Beta
        ##Adjust the Beta by weight vector 
        beta.ori2<-lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[c(2:(num.pc2+1)),]
        beta.adj2<-rbind(lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[1,],beta.ori2*weight2)
        ##Get the residual PCs
        yres2<-Y2-cbind(1,ypca2$x[,sign.pc2])%*%beta.adj2
        yres.pca2<-prcomp(yres2, retx = TRUE, center = T, scale = T, tol = NULL)
        Z2.est<-yres.pca2$x[,1:c(n.rpca2)]
      }
      ####
      TempZ<-Z2.est
      ####
    }
    ####update the number of hidden variables####
    n.rpca3<- num.sv(t(Y2),model.matrix(~1+X), B =20, method = "be", seed = NULL)
    Z2.est<-yres.pca2$x[,1:c(n.rpca3)]
    print(paste('Final number of hidden variables from weighted residual matrix',n.rpca3));
    ########
    out<-list(Z1.est=Z1.est,Z2.est=Z2.est,temp.y2=temp.y2,n.rpca3=n.rpca3)
  }
  return(out)
}





#########################################################################
########################################################################
#########################################################################
########################################################################


rPCA32_hidden<-function(Y,X,nComp,nIter,mCoverge,pLevelx,pLevelz){
  ###read in 
  Y<-as.matrix(Y);
  X<-as.matrix(X);
  n<-dim(Y)[1];m<-dim(Y)[2] 
  nx<-dim(X)[2]   ##number of covariates 
  n.rpca2<-nComp  ##number of top PCs taken as Z 
  n.iter<-nIter   ##repeated times of iterative step 
  pcutx<-pLevelx      ##significance level for x
  pcutz<-pLevelz      ##significance level for z
  pcutw<-0.05        ##significance level for weight vector 
  convcut<-mCoverge  ##converge cut for Y2
  if(n.rpca2==0){
    out<-NULL
  } 
  if(n.rpca2>=1){
    ################ Get Z1.est by All Features ################
    ####A.Detect significant PCs for observed covariates  
    ypca<-prcomp(Y, retx = TRUE, center = T, scale. = T, tol = NULL)
    pvalue.raw.x<-matrix(NA,(min(n,m)*0.5),nx)
    for ( ii in 1:(min(n,m)*0.5)){
      for( jj in 1:nx){
        pvalue.raw.x[ii,jj]<-summary(lm(X[,jj]~ypca$x[,ii]))$coef[2,4]
      }
    }
    pvalue.adj.x<-apply(pvalue.raw.x,2,function(xx){p.adjust(xx, method = "fdr")})
    sign.pc<-unique(which(pvalue.adj.x<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
    print(paste("On Y, selected PCs for X:"))
    print(sign.pc)
    num.pc<-length(sign.pc)
    if(num.pc==0){
      print('On Y, number of selected PC is 0, will take top PCs as Z')
      Z1.est<-ypca$x[,1:n.rpca2]   
      yres<-Y
      yres.pca<-ypca
    }
    if(num.pc!=0){
      ####B.Find the weight vector, suggested by Seunggeun
      ##Get p-value in model 1 and model 2; Get the weight based on the two p-values 
      pvalue.y.m1<-matrix(NA,num.pc,m)
      pvalue.y.m2<-matrix(NA,num.pc,m)
      weight<-matrix(NA,num.pc,m)
      for (kk in 1: num.pc){
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]],which=2)$pvalue
        pvalue.y.m1[kk,]<-p.adjust(temp, method = "fdr")
        ##
        temp<-F.mp(Y~ypca$x[,sign.pc[kk]]+X,which=2)$pvalue
        pvalue.y.m2[kk,]<-p.adjust(temp, method = "fdr")
        ##
        weight[kk,]<-0
        weight[kk,(pvalue.y.m1[kk,]<=pcutw)&(pvalue.y.m2[kk,]>=pcutw)]<-1
      }
      ####C. rPCA2 based on adjusted Beta
      ## Adjust beta by weight vector 
      beta.ori<-lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[c(2:(num.pc+1)),]
      beta.adj<-rbind(lm.mp(Y, formula=~ypca$x[,sign.pc], store.fitted = FALSE)$coef[1,],beta.ori*weight)
      ## Get the residual PCs
      yres<-Y-cbind(1,ypca$x[,sign.pc])%*%beta.adj
      yres.pca<-prcomp(yres, retx = TRUE, center = T, scale = T, tol = NULL)
      Z1.est<-yres.pca$x[,1:n.rpca2]
    }
    ################ Update Z1.est by Y2 (Subset of Features) ##############
    TempZ<-Z1.est
    ###Do iteration 
    for (gg in 1:n.iter){
      ####D. Get Y2: subset of Y with ZD>>XB
      pvalue.y2<-matrix(NA,2,m)
      temp<-F.mp(Y~X+TempZ,which=c(2:(1+nx)))$pvalue
      pvalue.y2[1,]<-p.adjust(temp, method = "fdr")
      temp<-F.mp(Y~TempZ+X,which=c(2:(1+n.rpca2)))$pvalue
      pvalue.y2[2,]<-p.adjust(temp, method = "fdr")
      temp<-(pvalue.y2[1,]>= pcutx)&(pvalue.y2[2,]<=pcutz)
      if(gg>1){
        diff.y2<-sum((temp.y2*1)!=(temp*1))
        print(paste("Difference (new Y2 and old Y2): ",diff.y2));
        if(diff.y2<=convcut){
          print('Y2 Converged, Done')
          break;
        }
      }
      temp.y2<-temp
      Y2<-Y[,temp]
      m2<-dim(Y2)[2]
      print(paste("Iter on Y2:",gg,", Dim :", m2,", n.hidden :", n.rpca2))
      ####E. Detect significant PCs for observed covariates on Y2
      ypca2<-prcomp(Y2, retx = TRUE, center = T, scale. = T, tol = NULL)
      pvalue.raw.x2<-matrix(NA,(min(n,m2)*0.5),nx)
      for ( ii in 1:(min(n,m2)*0.5)){
        for( jj in 1:nx){
          pvalue.raw.x2[ii,jj]<-summary(lm(X[,jj]~ypca2$x[,ii]))$coef[2,4]
        }
      }
      pvalue.adj.x2<-apply(pvalue.raw.x2,2,function(xx){p.adjust(xx, method = "fdr")})
      sign.pc2<-unique(which(pvalue.adj.x2<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
      print(paste('PCs from Y2 related to X: '));print(sign.pc2)
      num.pc2<-length(sign.pc2)
      ## If no selected PC, output the PCs
      if(num.pc2==0){
        print('On Y2, number of selected PC is 0, will take top PCs as Z')
        Z2.est<-ypca2$x[,1:n.rpca2]
        yres2<-Y2
        yres.pca2<-ypca2
      }
      ###If there is selected PC, do rPCA2 again  
      if(num.pc2>0){
        ## Get p-value in model 1 and model 2 
        ## Get the weight based on the two p-values 
        pvalue.y2.m1<-matrix(NA,num.pc2,m2)
        pvalue.y2.m2<-matrix(NA,num.pc2,m2)
        weight2<-matrix(NA,num.pc2,m2)
        for (kk in 1: num.pc2){
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]],which=2)$pvalue
          pvalue.y2.m1[kk,]<-p.adjust(temp, method = "fdr")
          ##
          temp<-F.mp(Y2~ypca2$x[,sign.pc2[kk]]+X,which=2)$pvalue
          pvalue.y2.m2[kk,]<-p.adjust(temp, method = "fdr")
          ##
          weight2[kk,]<-0
          weight2[kk,(pvalue.y2.m1[kk,]<=pcutw)&(pvalue.y2.m2[kk,]>=pcutw)]<-1
        }
        print(paste('Number of features removing any PC:',sum(apply(weight2,2,sum)>0)));
        ####F. rPCA2 based on adjusted Beta
        ##Adjust the Beta by weight vector 
        beta.ori2<-lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[c(2:(num.pc2+1)),]
        beta.adj2<-rbind(lm.mp(Y2, formula=~ypca2$x[, sign.pc2], store.fitted = FALSE)$coef[1,],beta.ori2*weight2)
        ##Get the residual PCs
        yres2<-Y2-cbind(1,ypca2$x[,sign.pc2])%*%beta.adj2
        yres.pca2<-prcomp(yres2, retx = TRUE, center = T, scale = T, tol = NULL)
        Z2.est<-yres.pca2$x[,1:c(n.rpca2)]
      }
      ####
      TempZ<-Z2.est
      ####
    }
    ####Remove effects of X on final Y2####
    for(hh in 1:n.iter){
      pvalue.raw.x.r<-matrix(NA,(min(n,m2)*0.5),nx)
      for ( ii in 1:(min(n,m2)*0.5)){
        for( jj in 1:nx){
          pvalue.raw.x.r[ii,jj]<-summary(lm(X[,jj]~yres.pca2$x[,ii]))$coef[2,4]
        }
      }
      pvalue.adj.x.r<-apply(pvalue.raw.x.r,2,function(xx){p.adjust(xx, method = "fdr")})
      sign.pc.r<-unique(which(pvalue.adj.x.r<=0.05,arr.ind=T)[,1]) #control FDR at 0.05 level
      print(paste("On R2, selected PCs for X:"))
      print(sign.pc.r)
      if(length(sign.pc.r)==0){break}
      ##update weight vector
      num.pc.r<-length(sign.pc.r)
      pvalue.y.m1<-matrix(NA,num.pc.r,m2)
      pvalue.y.m2<-matrix(NA,num.pc.r,m2)
      weight2<-matrix(NA,num.pc.r,m2)
      for (kk in 1: num.pc.r){
        temp<-F.mp(Y2~yres.pca2$x[,sign.pc.r[kk]],which=2)$pvalue
        pvalue.y.m1[kk,]<-p.adjust(temp, method = "fdr")
        ##
        temp<-F.mp(Y2~yres.pca2$x[,sign.pc.r[kk]]+X,which=2)$pvalue
        pvalue.y.m2[kk,]<-p.adjust(temp, method = "fdr")
        ##
        weight2[kk,]<-0
        weight2[kk,(pvalue.y.m1[kk,]<=pcutw)&(pvalue.y.m2[kk,]>=pcutw)]<-1
      }
      print(sum(apply(weight2,2,sum)>0))
      if(sum(apply(weight2,2,sum)>0)==0){break}
      ## update beta
      beta.ori2<-lm.mp(Y2, formula=~yres.pca2$x[,sign.pc.r], store.fitted = FALSE)$coef[c(2:(num.pc.r+1)),]
      beta.adj2<-rbind(lm.mp(Y2, formula=~yres.pca2$x[,sign.pc.r], store.fitted = FALSE)$coef[1,],beta.ori2*weight2)
      ## Update the residual PCs
      yres2<-Y2-cbind(1,yres.pca2$x[,sign.pc.r])%*%beta.adj2
      yres.pca2<-prcomp(yres2, retx = TRUE, center = T, scale = T, tol = NULL)
    }
    ####
    Z2.est<-yres.pca2$x[,1:c(n.rpca2)]
    Z3.est<-Z2.est
    for(hh in 1:1){
      beta.ori3<-lm.mp(Y2, formula=~X+Z3.est, store.fitted = FALSE)$coef[c(1:(nx+1)),]
      yres3<-Y2-cbind(1,X)%*%beta.ori3
      yres.pca3<-prcomp(yres3, retx = TRUE, center = T, scale = T, tol = NULL)
      Z3.est<-yres.pca3$x[,1:c(n.rpca2)]
    }
    ########
    out<-list(Z1.est=Z1.est,Z2.est=Z2.est,Z3.est=Z3.est,temp.y2=temp.y2)
  }
  return(out)
}
