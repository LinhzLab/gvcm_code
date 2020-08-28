###########################################
library(MASS)
linkFun <- function(Xbeta,case){
    sd=0.5;N.1 = length(Xbeta)
    if(case == 4){g = exp(Xbeta)/(1+exp(Xbeta));deriv.g = exp(Xbeta)/(1+exp(Xbeta))^2
    } else {g <- pnorm(Xbeta,0,sd);deriv.g <- dnorm(Xbeta,0,sd)}
   return(list(g,deriv.g))
}
genInput<-function(N,P)
{  
    X<-mvrnorm(N,rep(0,P),diag(P))
    U<-runif(N,0,1)
    Covariate<-list(x=X,u=U) 
  return(Covariate)
}
 #input_covariate <- genInput(N,P)
genData<-function(input_covariate,case){
  x<-input_covariate$x; u<-input_covariate$u
  P=ncol(x); N=nrow(x);y <- rep(0,N)
#######generate error########
  sd=0.5;error=rnorm(N,sd)
  beta.1=sin(pi*u);beta.2=cos(pi*u)
  bbeta=matrix(cbind(beta.1,beta.2),N,P)
  for(i in 1:N){bbeta[i,]=sign(bbeta[i,1])*bbeta[i,]}
  Xbeta<-rowSums(x*bbeta);
  if(case == 4){
    g=linkFun(Xbeta,case)[[1]];deriv.g=linkFun(Xbeta,case)[[2]];Variance=g*(1-g)
    for(i in 1:N){y[i]<-rbinom(1, size=1, g[i]);i = i+1};
  } else {for(i in 1:N){y[i]<-ifelse(Xbeta[i]+error[i]>=0,1,0);i = i+1};g=linkFun(Xbeta,case)[[1]];deriv.g=linkFun(Xbeta,case)[[2]];Variance=g*(1-g)
  }
  D<-list(Y=y,U = u,X=x,bbeta=bbeta,g=matrix(c(g,deriv.g),N,2),Variance=Variance,Xbeta=Xbeta)
  return(D)
}
#data <- genData(input_covariate,case)$g
