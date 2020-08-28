rm(list=ls())
library(splines)
library(splines2)
library(MASS)
library(ggplot2)
#N = 15000;P = 2;nk1 = 24;df = 25
N = 8000;P = 2;nk = 4;nk1 = 5;nk2= 4
for(s in 210:220){
  genInput<-function(N,P)
  {  
    X<-mvrnorm(N,rep(0,P),diag(P))
    U<-runif(N,0,1)
    Covariate<-list(x=X,u=U) 
    return(Covariate)
  }
  set.seed(s)
  input_covariate <- genInput(N,P)
  genData<-function(input_covariate)
  {#input<-INPUT[[sn]]
    x<-input_covariate$x; u<-input_covariate$u 
    P=ncol(x); N=nrow(x)
    #######generate error########
    beta.1=sin(0.5*pi*u); beta.2=cos(0.5*pi*u)
    bbeta=matrix(cbind(beta.1,beta.2),N,P)
    Xbeta<-rowSums(x*bbeta)
    g=5*pnorm(Xbeta,mean = 0, sd = 1); deriv.g=5*dnorm(Xbeta,mean = 0, sd = 1)
    Variance =(exp(-g)+1)^2
    error=rnorm(N,mean = 0,exp(-g)+1)
    y<-g+error
    #plot(g,Variance,ylim = c(0,0.4))
    #lines(g,(y-g)^2,type = 'p')
    plot(g,y,type = 'p')
    D<-list(Y=y,X=x,U=u,bbeta=bbeta,g=matrix(c(g,deriv.g),N,2),Variance=Variance,Xbeta=Xbeta)
    return(D)
  }
  data <- genData(input_covariate)
  #str(data)
  #plot(data$Xbeta,data$Y,type = 'p')
  Xbeta = data$Xbeta;y = data$Y;u = data$U;x=data$X;link = data$g[,1]
  #Bu <- bs(xbeta,df = df,intercept = TRUE)
  #inital_parameter <- rep(0.8,df)
  #parameter = optim(inital_parameter,objec_f)$par
  link_est <- function(y,alpha_2_0,X_Qta_T,V){
    B_star <- bSpline(X_Qta_T%*% alpha_2_0,df = nk1,intercept = TRUE)
    y_star <- y/sqrt(V);B_star_star <- apply(B_star,2,function(x){x/sqrt(V)})
    alpha_1 <- lm(y_star~B_star_star-1)$coefficients
    return(alpha_1)  
  }
  
  Beta_est <- function(y,alpha_2_0,alpha_1,X_Qta_T,V){
    B_star <- bSpline(X_Qta_T%*% alpha_2_0,df = nk1,intercept = TRUE)
    one_deriv <- deriv(B_star)%*%alpha_1
    covariate <- apply(X_Qta_T,2,function(x){x*one_deriv})
    intercept <- B_star %*% alpha_1- one_deriv*(X_Qta_T%*% alpha_2_0)
    y_star <- (y-intercept)/sqrt(V);covariate_star <- apply(covariate,2,function(x){x/sqrt(V)})
    alpha_2 <- lm(y_star~covariate_star-1)$coefficients
    return(alpha_2)
  }
  
  function_initial <- function(){
    N=length(y); P=ncol(x);bbeta = data$bbeta
    Bu <- bSpline(u,df = nk,intercept = TRUE);I_p = diag(1,2);k = 1
    #objective function, nk: the number of basis function, thus total parameter (p+1)*nk  
    # we use iterative algorithm to estimate bbeta and g
    # first we need to estimate initial value for bbeta
    X_Qta_T <- matrix(rep(0,N*P*nk),ncol = P*nk,byrow = T)
    for(i in 1:N){
      X_Qta_T[i,] <- x[i,]%*%kronecker(t(Bu[i,]),I_p)
    }
    #alpha_2_0 <- rep(0,2*nk)
    V <- rep(1,N)
    #V <-data$Variance
    alpha_2_0 <- lm(y~X_Qta_T-1)$coefficients
    #alpha_2_0[1:nk] <- lm(bbeta[,1]~Bu-1)$coefficients
    #alpha_2_0[-(1:nk)] <- lm(bbeta[,2]~Bu-1)$coefficients
    
    loss_0 <-10^5
    for(k in 1:70){
      alpha_1 <- link_est(y,alpha_2_0,X_Qta_T,V)
      alpha_2 <- Beta_est(y,alpha_2_0,alpha_1,X_Qta_T,V)
      ##add identify condition for alpha_2 bbeta(0) = 1
      scale <- sqrt((alpha_2[nk-3])^2+(alpha_2[nk-3+nk])^2)
      alpha_2[nk-3] <- alpha_2[nk-3] /scale; alpha_2[2*nk-3] <- alpha_2[2*nk-3] /scale
      #alpha_2[2*nk-3] <- 1;alpha_2[nk-3] <-0
      B_star <- bSpline(X_Qta_T%*% alpha_2,df = nk1,intercept = TRUE)
      loss_1 <- sum((y-B_star%*%alpha_1)^2)
      err = sum((alpha_2-alpha_2_0)^2)  
      if(err<10^-3|(loss_1-0.01-loss_0)>0){
        print("yes")
        break
      }
      #print(err)
      alpha_2_0 <- alpha_2
      loss_0 <- loss_1
    }
    return(c(alpha_1,alpha_2))
  }
  
  # Variance_initial <- function(y,nk2,mu){
  #   B_base <- bs(mu,df = nk2,intercept = TRUE)
  #   objec_f <- function(parameter){
  #     expression_0 <- y*y-B_base%*%parameter
  #     expression <- colSums(expression_0*expression_0)
  #   }
  #   in_parameter <- rep(0.1,nk2)
  #   parameter <- optim(in_parameter,objec_f)$par
  #   variance <-B_base%*%parameter-mu^2 
  #   return(variance)
  # }
  
  Variance_initial <- function(y,nk2,mu){
    Y <-(y-mu)^2;
    B_base <- bs(mu,df = nk2,intercept = TRUE)
    theta <- lm(Y~B_base-1)$coefficients
    variance <-B_base%*%theta
    return(variance)
  }
  
  parameter_0 <- function_initial()
  alpha_2 <- parameter_0[-(1:nk1)]
  dim(alpha_2) <- c(nk,P)
  Bu <- bs(u,df =nk ,intercept = TRUE)
  bbeta_0 <- Bu%*%alpha_2
  #Ident=sqrt(sum(bbeta_0[which.min(abs(u-mean(u))),1:P]^2));for(i in 1:N){bbeta_0[i,]=sign(bbeta_0[i,1])*bbeta_0[i,]/Ident}
  xbeta_0 = rowSums(x*bbeta_0)
  Bx = bs(xbeta_0,df = nk1,intercept = TRUE)
  f_hat <- Bx%*%parameter_0[1:nk1]
  vari_initial <- Variance_initial(y,nk2,f_hat)
  #vari_initial<- rep(1,N)
  print(mean(sapply(data$bbeta[,1]-bbeta_0[,1],function(x){x^2})))
  #mean(sapply(data$bbeta[,2]-y_beta_1,function(x){x^2}))
  print(mean(sapply(data$bbeta[,2]-bbeta_0[,2],function(x){x^2})))
  print(mean(sapply(data$g[,1]-f_hat,function(x){x^2})))#print(mean(sapply(data$Variance-vari_initial,function(x){x^2})))
  print(mean(sapply(data$Variance-vari_initial,function(x){x^2})))
  print(s)
  # plot(Xbeta[order(Xbeta)],link[order(link)],xlim = c(-5,5),ylim = c(0,1.2),type = 'l',cex=1.5,font=1)
  # lines(xbeta_0,f_hat,ylim = c(0,1.2),type = 'p',cex=0.1,font=1)
  #plot(data$U,bbeta_0[,1],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  #lines(data$U,data$bbeta[,1],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  #plot(data$U,bbeta_0[,2],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  #lines(data$U,data$bbeta[,2],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  
  N=length(y); P=ncol(x);bbeta = data$bbeta
  Bu <- bSpline(u,df = nk,intercept = TRUE);I_p = diag(1,2);k = 1
  #objective function, nk: the number of basis function, thus total parameter (p+1)*nk  
  # we use iterative algorithm to estimate bbeta and g
  # first we need to estimate initial value for bbeta
  V <- vari_initial
  alpha_2_0 <- parameter_0[-(1:nk1)]
  #alpha_2_0[1:nk] <- lm(bbeta[,1]~Bu-1)$coefficients
  #alpha_2_0[-(1:nk)] <- lm(bbeta[,2]~Bu-1)$coefficients
  loss_0 <-10^5
  X_Qta_T <- matrix(rep(0,N*P*nk),ncol = P*nk,byrow = T)
  for(i in 1:N){
    X_Qta_T[i,] <- x[i,]%*%kronecker(t(Bu[i,]),I_p)
  }
  for(k in 1:70){
    alpha_1 <- link_est(y,alpha_2_0,X_Qta_T,V)
    alpha_2 <- Beta_est(y,alpha_2_0,alpha_1,X_Qta_T,V)
    ##add identify condition for alpha_2 bbeta(0) = 1
    scale <- sqrt((alpha_2[nk-3])^2+(alpha_2[nk-3+nk])^2)
    alpha_2[nk-3] <- alpha_2[nk-3] /scale; alpha_2[2*nk-3] <- alpha_2[2*nk-3] /scale
    #alpha_2[2*nk-3] <- 1;alpha_2[nk-3] <-0
    B_star <- bSpline(X_Qta_T%*% alpha_2,df = nk1,intercept = TRUE)
    loss_1 <- sum((y-B_star%*%alpha_1)^2)
    err = sum((alpha_2-alpha_2_0)^2)  
    if(err<10^-3|(loss_1-0.01-loss_0)>0){
      print("yes")
      break
    }
    #print(err)
    alpha_2_0 <- alpha_2
    loss_0 <- loss_1
  }
  dim(alpha_2) <- c(nk,P)
  Bu <- bs(u,df =nk ,intercept = TRUE)
  bbeta_0 <- Bu%*%alpha_2
  #Ident=sqrt(sum(bbeta_0[which.min(abs(u-mean(u))),1:P]^2));for(i in 1:N){bbeta_0[i,]=sign(bbeta_0[i,1])*bbeta_0[i,]/Ident}
  xbeta_0 = rowSums(x*bbeta_0)
  Bx = bs(xbeta_0,df = nk1,intercept = TRUE)
  f_hat <- Bx%*% alpha_1
  vari_initial <- Variance_initial(y,nk2,f_hat)
  
  
  print(mean(sapply(data$bbeta[,1]-bbeta_0[,1],function(x){x^2})))
  #mean(sapply(data$bbeta[,2]-y_beta_1,function(x){x^2}))
  print(mean(sapply(data$bbeta[,2]-bbeta_0[,2],function(x){x^2})))
  print(mean(sapply(data$g[,1]-f_hat,function(x){x^2})))#print(mean(sapply(data$Variance-vari_initial,function(x){x^2})))
  print(mean(sapply(data$Variance-vari_initial,function(x){x^2})))
  plot(f_hat,vari_initial,xlim = c(0,1),ylim = c(0,0.4),color = 'blue',type = 'p',cex=0.5,font=1)
  lines(f_hat,data$Variance ,xlim = c(0,1),ylim = c(0,0.4),type = 'p',cex=0.5,font=1)
  print(s)
  # plot(Xbeta[order(Xbeta)],link[order(link)],xlim = c(-5,5),ylim = c(0,1.2),type = 'l',cex=1.5,font=1)
  # lines(xbeta_0,f_hat,ylim = c(0,1.2),type = 'p',cex=0.1,font=1)
  plot(data$U,bbeta_0[,1],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  lines(data$U,data$bbeta[,1],xlim = c(0,1),ylim = c(-2,2),type = 'p',cex=0.5,font=1)
  #plot(f_hat,vari_initial,color = 'blue',type = 'p',cex=0.5,font=1)
  #lines(f_hat,(exp(-f_hat)+2.5)^2,xlim = c(0,1),ylim = c(0,5),type = 'p',cex=0.5,font=1)
  
}
#parameter <- optim(parameter_0,Quasi_initial)$par

#alpha_2 <- parameter[-(1:nk1)]
#dim(alpha_2) <- c(df,P)
#Bu <- bs(u,df = df,intercept = TRUE)
#bbeta <- Bu%*%alpha_2
##Ident=sqrt(sum(bbeta[which.min(abs(u-mean(u))),1:P]^2));for(i in 1:N){bbeta[i,]=sign(bbeta[i,1])*bbeta[i,]/Ident}
#xbeta = rowSums(x*bbeta)
#Bx = bs(xbeta,df = nk1,intercept = TRUE)
#f_hat <- Bx%*%parameter[1:nk1]
#y_beta = bbeta[,1]
#y_beta_1 = bbeta[,2]
#vari_est <- Variance_initial(y,nk2,f_hat)


