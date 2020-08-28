###########################################
library(MASS)
linkFun <- function(Xbeta,case){
    if(case == 1){g = Xbeta;deriv.g = rep(1,N)
    } else if (case == 2){g <- Xbeta^2;deriv.g <- 2*Xbeta
    } else if (case == 3){g <- sin(2*Xbeta);deriv.g <- 2*cos(2*Xbeta) 
    } else {g <- pnorm(Xbeta,mean = 0, sd = 1);g <- dnorm(Xbeta,mean = 0, sd = 1)}
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
genData<-function(input_covariate)
{#input<-INPUT[[sn]]
  x<-input_covariate$x; u<-input_covariate$u 
  P=ncol(x); N=nrow(x)
#######generate error########
  sd=0.1;error=rnorm(N,mean = 0,sd)
  beta.1=u^2+1; beta.2=cos(pi*u)*cos(pi*u)+rep(0.5,N);beta.3=2*sin(pi*u)*sin(pi*u)-rep(0.5,N)
  bbeta=matrix(cbind(beta.1,beta.2,beta.3),N,P)
for(i in 1:N){bbeta[i,]=sign(bbeta[i,1])*bbeta[i,]/sqrt(sum(bbeta[i,1:P]^2))}
  Xbeta<-rowSums(x*bbeta);g=linkFun(Xbeta,case)[[1]]; deriv.g=linkFun(Xbeta,case)[[2]];Variance=rep(sd^2,N); 
  y<-g+error
  D<-list(Y=y,X=x,U=u,bbeta=bbeta,g=matrix(c(g,deriv.g),N,2),Variance=Variance,error=error,Xbeta=Xbeta)
  return(D)
}
#data <- genData(input_covariate)$g
