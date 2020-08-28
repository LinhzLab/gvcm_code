###########################################
library(MASS)
linkFun <- function(Xbeta){
  g <- 5*pnorm(Xbeta,mean = 0, sd = 1);g.deriv <- 5*dnorm(Xbeta,mean = 0, sd = 1)
  return(list(g,g.deriv))
}
genInput<-function(N,P)
{  
  X<-mvrnorm(N,rep(0,P),diag(P))
  U<-runif(N,0,1)
  Covariate<-list(x=X,u=U) 
  return(Covariate)
}

genData<-function(input_covariate)
{#input<-INPUT[[sn]]
  x<-input_covariate$x; u<-input_covariate$u 
  P=ncol(x); N=nrow(x)
  #######generate error########
  beta.1=sin(0.5*pi*u); beta.2=cos(0.5*pi*u)
  bbeta=matrix(cbind(beta.1,beta.2),N,P)
  for(i in 1:N){bbeta[i,1]=sign(bbeta[i,1])*bbeta[i,1]}
  bbeta[which.min(u),] <- bbeta[which.min(u),]/sqrt(sum(bbeta[which.min(u),]*bbeta[which.min(u),]))
  Xbeta<-rowSums(x*bbeta)
  g=linkFun(Xbeta)[[1]]; deriv.g=linkFun(Xbeta)[[2]];
  Variance =(exp(-g)+1)^2
  error=rnorm(N,mean = 0,exp(-g)+1)
  y<-g+error
  D<-list(Y=y,X=x,U=u,bbeta=bbeta,g=matrix(c(g,deriv.g),N,2),Variance=Variance,error=error,Xbeta=Xbeta)
  return(D)
}

