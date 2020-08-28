  library(mgcv)
  library(fda)

  CV1=function(N,KF){
  z=rep(1:KF,ceiling(N/KF))[1:N]
  set.seed(44);z=sample(z,N)
  mm=list();for (i in 1:KF) mm[[i]]=(1:N)[z==i]
  return(mm)}

    P_prob <- function(x){
   if(x<=710){pre_p = exp(x)/(1+exp(x))
   }else{ pre_p = 1}
    return(pre_p)
  }
   P_prob_deriv <- function(x){
   if(x<=710){pre_p = exp(x)/(1+exp(x))^2
   }else{ pre_p = 0}
    return(pre_p)
  } 

  K<-function(h,u,mu)
  { N=length(u)
    Knel<-1/h*3/4*(1-((u-mu)/h)^2)*ifelse(abs((u-mu)/h)<=1,1,0)
    return(Knel)
  }
# estimate \beta function vc:beta
 vc.hf<-function(h,Data,Beta,dp)
   { #h = c(0.15,0.2);Data = DATA[[141]];Beta = Para.initial;Beta_dp = BETA0_dp;dp = U.fix
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    vc0<-cbind(Beta$bbeta,matrix(0,N,P)); g<-Beta$g[,1]; dg<-Beta$g[,2]; V<-Beta$V
    vc1=matrix(NA,length(dp),2*P)
    for(i in 1:length(dp))
    {
      M.U<-c(u-dp[i]); GMU<-matrix(c(x,M.U*x),N,2*P)
      #W1=c(K(h,u,dp[i])*dg/V) ; W2=c(K(h,u,dp[i])*dg^2/V) 
      #K_v <-cbind(K(h[1],u,dp[i]),K(h[1],u,dp[i]),K(h[2],u,dp[i]))
      #K_v <- cbind(K_v,K_v)
      W1=c(K(h,u,dp[i])*dg/(V+10^-5)) ; W2=c(K(h,u,dp[i])*dg^2/(V+10^-5)) 
      yita=y-g+dg*c(rowSums(x*vc0[,1:P]))
      PJ=t(GMU)%*%(W2*GMU) 
      vc1[i,]=solve(PJ+diag(10^-5,2*P))%*%t(GMU)%*%(W1*yita)#} 
      #yita=y-g+dg*c(rowSums(x*vc0[,1:P]))
      #PJ=t(GMU)%*%(GMU*(K_v*c(dg^2/(V+10^-5))))
      #vc1[i,]=solve(PJ+diag(10^-5,2*P))%*%t(GMU*K_v)%*%(c(dg/(V+10^-5))*yita)#} 
    }
    return(vc1)
  }  
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

#derive g and \dot g function by (2.6) in paper 
  link.kernel<-function(h,data,beta_hat,Vari,dp){
    y<-data$Y; x <- data$X; N=length(y) 
    linkh=matrix(NA,length(dp),2)
    for(i in 1:length(dp)){
      W=c(K(h,rowSums(x*beta_hat),dp[i])/Vari);z<-cbind(rep(1,N),rowSums(x*beta_hat)-dp[i])
      PJ=t(z)%*%(z*W)
      linkh[i,]=solve(PJ+diag(10^-4,2))%*%t(z)%*%(W*y)#}
      i = i+1
      }
      return(linkh)
  }
  

v.hf<-function(h,Data,link.est,link.fix){#h = h[3];Data<-DATA[[sn]];g = gh[,1];dp<-g.fix
  y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
  yita=(y-link.est)^2
  vh=matrix(NA,length(link.fix),2)
  for(i in 1:length(link.fix))
  {
    W=c(K(h,link.est,link.fix[i])); z<-cbind(rep(1,N),link.est-link.fix[i])
    PJ3=t(z)%*%(W*z); vh[i,]=solve(PJ3+diag(10^-4,2))%*%t(z)%*%(W*yita)
  }
  return(vh)
}
  
spline.initial<-function(sn){  
   y<-DATA[[sn]]$Y; x <- DATA[[sn]]$X; u<- DATA[[sn]]$U; N=length(y); P=ncol(x);Xbeta = DATA[[sn]]$Xbeta
   Bu <- bSpline(u,df = nk,intercept = TRUE);I_p = diag(1,3);i = 1;k = 1
    #objective function, nk: the number of basis function, thus total parameter (p+1)*nk  
    # we use iterative algorithm to estimate bbeta and g
    # first we need to estimate initial value for bbeta
   X_Qta_T <- matrix(rep(0,N*P*nk),ncol = P*nk,byrow = T)
   for(i in 1:N){
     X_Qta_T[i,] <- x[i,]%*%kronecker(I_p,t(Bu[i,]))
   }
   #alpha_2_0 <- rep(0,2*nk)
   V <- rep(0.01,N)
   #V <-data$Variance
   alpha_2_0 <- lm(y~X_Qta_T-1)$coefficients
   #alpha_2_0[1:nk] <- lm(bbeta[,1]~Bu-1)$coefficients
   #alpha_2_0[-(1:nk)] <- lm(bbeta[,2]~Bu-1)$coefficients
   loss_0 <-10^5
   for(k in 1:70){
     alpha_1 <- link_est(y,alpha_2_0,X_Qta_T,V)
     alpha_2 <- Beta_est(y,alpha_2_0,alpha_1,X_Qta_T,V)
     ##add identify condition for alpha_2 bbeta(0) = 1
     #scale <- sqrt((alpha_2[1])^2+(alpha_2[1+nk])^2+(alpha_2[1+2*nk])^2)
     #alpha_2[1] <- alpha_2[1] /scale; alpha_2[1+nk] <- alpha_2[1+nk] /scale
     #alpha_2[1+2*nk] <- alpha_2[1+2*nk] /scale
     B_star <- bSpline(X_Qta_T%*% alpha_2,df = nk1,intercept = TRUE)
     loss_1 <- sum((y-B_star%*%alpha_1)^2)
     err = sum((alpha_2-alpha_2_0)^2)  
     if(err<10^-3|(loss_1-0.1-loss_0)>0){
       print("yes")
       break
     }
     #print(err)
     alpha_2_0 <- alpha_2
     loss_0 <- loss_1
   }
  beta_hat <- Bu%*%matrix(alpha_2,nrow = nk)
  for(i in 1:N){beta_hat[i,]=sign(beta_hat[i,1])*beta_hat[i,]/sqrt(sum(beta_hat[i,1:P]^2))}
   xbeta<-rowSums(x*beta_hat); Bx<-bs(xbeta,df = nk1,intercept = TRUE)
   g_hat<-Bx%*%alpha_1
   dp=seq(range(xbeta)[1],range(xbeta)[2],diff(c(range(xbeta)))/Points)
   g.deriv <- link.kernel(h[2],DATA[[sn]],beta_hat,V,dp)[,2]
   g.deriv = approx(dp,g.deriv,rule = 2,xout=xbeta)$y
   #beta_hat[1:20,]
   #DATA[[sn]]$bbeta[1:20,]
   sp.est = data.frame(beta = beta_hat, m =as.vector(g_hat), dm = g.deriv, Vari = V)
   #print(paste0(mean(sapply(data$bbeta[,1]-sp.est[,1],function(x){x^2})), mean(sapply(data$bbeta[,2]-sp.est[,2],function(x){x^2})), mean(sapply(data$bbeta[,3]-sp.est[,3],function(x){x^2})),sep='_:_',collapse = "; "))
   return(sp.est)
}

calculate_gvcm_parallel<-function(sn){
    x<-DATA[[sn]]$X;U<-DATA[[sn]]$U;P=ncol(DATA[[sn]]$X);N=length(DATA[[sn]]$U)
    Para.initial<-list(bbeta=gm1[[sn]][,1:P],g=cbind(gm1[[sn]][,P+1],g.deriv = gm1[[sn]][,P+2]),V = gm1[[sn]][,P+3])
    U.fix<-seq(range(U)[1],range(U)[2],diff(range(U))/Points) 
    bbeta.fix.est<-vc.hf(h[1],DATA[[sn]],Para.initial,U.fix)   
    bbeta.fix.1.est<-matrix(NA,(Points+1),P);bbeta.est<-matrix(NA,N,2*P)
    bbeta_true_fix <- matrix(NA,(Points+1),P)  
    for(k in 1:P)
    {
      bbeta.est[,k]<-approx(U.fix,bbeta.fix.est[,k],xout=U,rule=2)$y;bbeta.est[,P+k]<-approx(U.fix,bbeta.fix.est[,P+k],xout=U,rule=2)$y
     } 
    for(i in 1:N){bbeta.est[i,]=sign(bbeta.est[i,1])*bbeta.est[i,]/sqrt(sum(bbeta.est[i,1:P]^2))}# vc1(0)=1 
     U.fix.1<-seq(quantile(U,probs = 0.1),quantile(U,probs = 0.9),diff(c(quantile(U,probs = 0.1),quantile(U,probs = 0.9)))/Points)
     for(k in 1:P)
     {
       bbeta.fix.1.est[,k] <- approx(U,bbeta.est[,k],xout=U.fix.1,rule=2)$y
     }
       for(k in 1:P)
     {
       bbeta_true_fix[,k] <- approx(U,DATA[[sn]]$bbeta[,k],xout=U.fix.1,rule=2)$y
     }
    Xbeta.est=rowSums(bbeta.est[,1:P]*x)
    #采用quantile来取点，K估计有效，不会周围没有点
    Xbeta.fix<-seq(quantile(Xbeta.est,probs = 0.1),quantile(Xbeta.est,probs = 0.9),diff(c(quantile(Xbeta.est,probs = 0.1),quantile(Xbeta.est,probs = 0.9)))/(Points+100))
    g_true_fix <- linkFun(Xbeta.fix,case)[[1]]
    g.fix.est<-link.kernel(h[2],DATA[[sn]],bbeta.est[,1:P],Para.initial$V,Xbeta.fix);link.fix.est <- g.fix.est[,1]
    link.est<-approx(Xbeta.fix,g.fix.est[,1],xout=DATA[[sn]]$Xbeta,rule=2)$y 
    bbeta.est.P<-matrix(NA,(Points+1),P)
    for(k in 1:P)
    {#bbeta.est at Points from dp1 in order to plot
       bbeta.est.P[,k]<-approx(U.fix,bbeta.fix.est[,k],xout=dp1,rule=2)$y
    }
    link.est.P<-matrix(NA,Points+1,1)
    link.est.P<-approx(Xbeta.fix,link.fix.est,xout=dp2,rule=2)$y
    link.fix<-seq(quantile(link.est,probs = 0.1),quantile(link.est,probs = 0.9),diff(c(quantile(link.est,probs = 0.1),quantile(link.est,probs = 0.9)))/Points)
    Variance.fix.est <-try(v.hf(h[3],DATA[[sn]],link.est,link.fix)[,1],silent = TRUE)
    if('try-error' %in% class(Variance.fix.est)){
      next
    }else{
      if(!is.null(Variance.fix.est[Variance.fix.est<0])){
        Variance.fix.est[Variance.fix.est<0] <- min(Variance.fix.est[Variance.fix.est>0])
      }
    Variance.est.P<-approx(link.fix,Variance.fix.est,xout=dp3,rule=2)$y 
    par_bbeta_var <- cbind(bbeta.fix.1.est,Variance.fix.est,bbeta_true_fix,Variance_true =rep(0.01,(Points+1)))
    par_g <- cbind(link.fix.est,g_true_fix)
    parameter_fix <- cbind(bbeta.est.P,link.est.P,Variance.est.P) 
    }
    return(list(result_g = par_g,result_bbeta_var = par_bbeta_var,result.fix = parameter_fix))
}
##############################
# spline.initial<-function(sn){  
#   y<-DATA[[sn]]$Y; x <- DATA[[sn]]$X; u<- DATA[[sn]]$U; N=length(y); P=ncol(x)  
#   Bu <- bSpline(u,df = nk,intercept = TRUE) 
#   #objective function, nk: the number of basis function, thus total parameter (p+1)*nk    
#   objec_f <- function(parameter){
#     alpha_2 <- parameter[-(1:nk1)]
#     dim(alpha_2) <- c(nk,P)
#     xu <- rowSums(x*(Bu%*%alpha_2))
#     Bxu_1 <- bSpline(xu,df = nk1,intercept = TRUE)
#     expression_1 <- y-Bxu_1%*%parameter[1:nk1]
#     expression <- colSums(expression_1*expression_1)
#   }
#   inital_parameter <- rep(0.7,(P*nk+nk1))
#   parameter <- optim(inital_parameter,objec_f)$par
#   beta_hat <- Bu%*%t(matrix(parameter[-(1:nk1)],nrow = P,byrow = T))
#   Ident=sqrt(sum(beta_hat[which.min(abs(u-mean(u))),1:P]^2));for(i in 1:N){beta_hat[i,]=sign(beta_hat[i,1])*beta_hat[i,]/Ident}
#   beta.initial <- beta_hat
#   alpha.1.initial <- parameter[1:nk1]
#   alpha.2.initial <- parameter[-(1:nk1)]
#   tol = 10^2;k = 1
#   while(tol>0.1&&k<=50){
#     xu <- rowSums(x*beta.initial)
#     Bxu_1 <- bSpline(xu,df = nk1,intercept = TRUE)
#     alpha.1.est <- solve(t(Bxu_1)%*%Bxu_1+diag(10^-5,nk1))%*%t(Bxu_1)%*%y
#     alpha.1.initial <- alpha.1.est
#     #g <- Bxu_1%*%alpha.1.est
#     #Data$g[,1]-g
#     objec_alpha.2 <- function(alpha.2){
#       #     alpha.2 <- alpha.2.initial
#       dim(alpha.2) <- c(nk,P) #alpha.2;alpha.2.initial
#       xu <- rowSums(x*(Bu%*%alpha.2))
#       Bxu_1 <- bSpline(xu,df = nk1,intercept = TRUE)
#       expression_1 <- y-Bxu_1%*%alpha.1.initial
#       expression <- colSums(expression_1*expression_1)
#     }
#     alpha.2.initial <- optim(alpha.2.initial,objec_alpha.2)$par  
#     beta_hat <- Bu%*%matrix(alpha.2.initial,ncol = P,byrow = T)
#     Ident=sqrt(sum(beta_hat[which.min(abs(u-mean(u))),1:P]^2));for(i in 1:N){beta_hat[i,1]=sign(beta_hat[i,1])*beta_hat[i,1]/Ident}
#     tol = sum(abs(beta_hat-beta.initial))
#     print(tol)
#     beta.initial <- beta_hat
#     k = k+1
#   }
#   print(paste("the last tol is ",tol,seq=""))
#   beta_hat; xu<-rowSums(x*beta_hat); Bxu_1<-bs(xu,df = nk1,intercept = TRUE)
#   g_hat<-Bxu_1%*%alpha.1.est;parameter<-c(alpha.1.initial,alpha.2.initial)
#   return(list(beta_hat,g_hat,parameter))
#  } 