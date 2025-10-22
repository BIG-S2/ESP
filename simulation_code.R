
#simulation use only the real data dimension
library(lfmm)
library(pracma)
lfmm<-function(y,x,k,lambda=10^-8)
{
  x=as.matrix(x)
  aaa1=lfmm:::compute_eigen_svd(x)
  y_rota=t(aaa1$Q)%*%as.matrix(y)
  y_rota[1:dim(x)[2],]=y_rota[1:dim(x)[2],]*lambda
  space_rota=as.matrix(svd(y_rota)$u[,1:k])
  space_rota[1:dim(x)[2],]=1/lambda*space_rota[1:dim(x)[2],]
  conf_space=aaa1$Q%*%space_rota
  return(conf_space)
}

lfmm_getw<-function(y,x,k,lambda=10^-8)
{
  x=as.matrix(x)
  aaa1=lfmm:::compute_eigen_svd(x)
  y_rota=t(aaa1$Q)%*%as.matrix(y)
  y_rota[1:dim(x)[2],]=y_rota[1:dim(x)[2],]*lambda
  tt=svd(y_rota)
  space_rota=as.matrix(tt$u[,1:k]%*%diag(tt$d[1:k])%*%t(tt$v[,1:k]))
  space_rota[1:dim(x)[2],]=1/lambda*space_rota[1:dim(x)[2],]
  w=aaa1$Q%*%space_rota
  return(w)
}
