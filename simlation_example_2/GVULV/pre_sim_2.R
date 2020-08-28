  library(mgcv)
  library(fda)

  CV1=function(N,KF){
  z=rep(1:KF,ceiling(N/KF))[1:N]
  set.seed(44);z=sample(z,N)
  mm=list();for (i in 1:KF) mm[[i]]=(1:N)[z==i]
  return(mm)}
  K<-function(h,u,mu)
  { N=length(u)
    Knel<-1/h*3/4*(1-((u-mu)/h)^2)*ifelse(abs((u-mu)/h)<=1,1,0)
    return(Knel)
  }
# estimate \beta function vc:beta
vc.hf<-function(h,Data,Beta,dp)
  { #h = 1;Data = DATA[[sn]];Beta = Para.initial;dp = U.fix
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    vc0<-cbind(Beta$bbeta,matrix(0,N,P)); g<-Beta$g[,1]; dg<-Beta$g[,2]; V<-Beta$V
    vc1=matrix(NA,length(dp),2*P)
    for(i in 1:length(dp))
    {
      M.U<-c(u-dp[i]); GMU<-matrix(c(x,M.U*x),N,2*P)
      W1=c(K(h,u,dp[i])*dg/(V+10^-5)) ; W2=c(K(h,u,dp[i])*dg^2/(V+10^-5)) 
      yita=y-g+dg*c(rowSums(x*vc0[,1:P]))
      PJ=t(GMU)%*%(W2*GMU) 
      vc1[i,]=solve(PJ+diag(10^-4,2*P))%*%t(GMU)%*%(W1*yita)#} 
    }
    return(vc1)
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
#derive variance 
v.hf<-function(h,Data,link.est,link.fix)
  {#h = h[3];Data<-DATA[[sn]];g = gh[,1];dp<-g.fix
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    yita=(y)^2
    vh=matrix(NA,length(link.fix),2)
    for(i in 1:length(link.fix))
    {
      W=c(K(h,link.est,link.fix[i])); z<-cbind(rep(1,N),link.est-link.fix[i])
      PJ3=t(z)%*%(W*z); vh[i,]=solve(PJ3+diag(10^-4,2))%*%t(z)%*%(W*yita)
    }
    return(vh)
}
Variance_initial <- function(data,nk2,mu){
  y<-data$Y	
  B_base <- bs(mu,df = nk2,intercept = TRUE)
  objec_f <- function(parameter){
    expression_0 <- y*y-B_base%*%parameter
    expression <- colSums(expression_0*expression_0)
  }
#########choose inital paramter by true model ######
  True_V <- function(parameter){
    QtaV <- data$Variance+(data$g[,1])*(data$g[,1])
    B_base <- bs(data$g[,1],df = nk2,intercept = TRUE)
    expression_0 <- QtaV-B_base%*%parameter
    expression <- colSums(expression_0*expression_0)
  }
  in_parameter <- optim(rep(0,nk2),True_V)$par
  parameter <- optim(in_parameter,objec_f)$par
  variance <-B_base%*%parameter-mu^2 
  return(variance)
}
#use spline method to get the inital parameters of \beta, \dot g and  g; P: dimension of X
spline.initial<-function(sn){  
  y<-DATA[[sn]]$Y; x <- DATA[[sn]]$X; u<- DATA[[sn]]$U; N=length(y); P=ncol(x);Xbeta = DATA[[sn]]$Xbeta
  Bu <- bs(u,df = nk,intercept = TRUE) 
#objective function, nk: the number of basis function, thus total parameter (p+1)*nk    
  function_initial <- function(parameter){
    alpha_2 <- parameter[-(1:nk1)]
    dim(alpha_2) <- c(nk,P)
    Bu <- bs(u,df = nk,intercept = TRUE)
    xbeta = rowSums(x*(Bu%*%alpha_2))
    Bx = bs(xbeta,df = nk1,intercept = TRUE)
    expression_1 <- y-Bx%*%parameter[1:nk1]
    expression = sum(expression_1*expression_1)
  }
  Beta_initial <- function(parameter){
    dim(parameter) <- c(nk,P)
    Bu <- bs(u,df = nk,intercept = TRUE)
    expression_1 <- Xbeta-rowSums(x*(Bu%*%parameter))
    expression = sum(expression_1*expression_1)
  }
  link_initial <- function(parameter){
    Bx <- bs(Xbeta,df = nk1,intercept = TRUE)
    expression_1 <- y-Bx%*%parameter
    expression = sum(expression_1*expression_1)
  }
  alpha.2.inital <- optim(c(rep(0.1,P*nk)),Beta_initial)$par
  alpha.1.inital <- optim(c(rep(0.1,nk1)),link_initial )$par
  parameter_0 <- optim(c(alpha.1.inital,alpha.2.inital),function_initial)$par
  beta_hat <- Bu%*%t(matrix(parameter_0[-(1:nk1)],nrow = P,byrow = T))
  xbeta<-rowSums(x*beta_hat); Bx<-bs(xbeta,df = nk1,intercept = TRUE)
  g_hat<-Bx%*%parameter_0[1:nk1]
  vari_initial <-g_hat*(1-g_hat) 
  vari_initial[vari_initial<0] <- min(vari_initial[vari_initial>0])
  dp=seq(range(xbeta)[1],range(xbeta)[2],diff(c(range(xbeta)))/Points)
  g.deriv <- link.kernel(h[2],DATA[[sn]],beta_hat,rep(1,N),dp)[,2]
  g.deriv = approx(dp,g.deriv,rule = 2,xout=DATA[[sn]]$Xbeta)$y
  #plot_data <- cbind(xbeta = xbeta,g = g_hat);plot_data <- plot_data[order(plot_data[,1]),]
  #plot(plot_data[,1],plot_data[,2]);dev.off()
  sp.est = data.frame(beta = beta_hat, m =as.vector(g_hat), dm = g.deriv, Vari = as.vector(vari_initial))
  #print(paste0(mean(abs(sp.est[,1]-DATA[[1]]$bbeta[,1])), sep=',', mean(abs(sp.est[,2]-DATA[[1]]$bbeta[,2])),sep=',', mean(abs(sp.est[,3]-DATA[[1]]$g[,1])), sep=',', mean(abs(sp.est[,5]-DATA[[1]]$Variance)),sep=',', mean(abs(sp.est[,4]-DATA[[1]]$g[,2]),na.rm = TRUE),collapse = "; "))
  return(sp.est)
}
#gm1 = list()
#for (i in 5:10){
#	    gm1[[i]] = spline.initial(i);i = i+1}

calculate_gvcm_parallel<-function(sn){
  x<-DATA[[sn]]$X;U<-DATA[[sn]]$U;P=ncol(DATA[[sn]]$X);N=length(DATA[[sn]]$U);Xbeta = DATA[[sn]]$Xbeta
  Para.initial<-list(bbeta=gm1[[sn]][,1:P],g=cbind(gm1[[sn]][,P+1],g.deriv = gm1[[sn]][,P+2]),V = gm1[[sn]][,P+3])
####采用quantile来取点，K估计有效，不会周围没有点
####Fix point to estimate bbeta(U)####
  U.fix<-seq(quantile(U,probs = 0.01),quantile(U,probs = 0.99),diff(c(quantile(U,probs = 0.01),quantile(U,probs = 0.99)))/Points)    
  bbeta.fix.1.est<-matrix(NA,(Points+1),P)
  bbeta.est<-matrix(NA,N,2*P);bbeta_true_fix <- matrix(NA,(Points+1),P)     
  bbeta.fix.est<-vc.hf(h[1],DATA[[sn]],Para.initial,U.fix)       
  for(k in 1:P){
      bbeta.est[,k]<-approx(U.fix,bbeta.fix.est[,k],xout=U,rule=2)$y;bbeta.est[,P+k]<-approx(U.fix,bbeta.fix.est[,P+k],xout=U,rule=2)$y
  } 
####Fix point to evaluate the estimation performance####
  print(paste0(mean(abs(bbeta.est[,1]-DATA[[1]]$bbeta[,1])),sep=',',mean(abs(bbeta.est[,2]-DATA[[2]]$bbeta[,2]))))  

  U.fix.1<-seq(quantile(U,probs = 0.1),quantile(U,probs = 0.9),diff(c(quantile(U,probs = 0.1),quantile(U,probs = 0.9)))/Points)
  for(k in 1:P)
  {
    bbeta.fix.1.est[,k] <- approx(U,bbeta.est[,k],xout=U.fix.1,rule=2)$y
  }
#####Find the true function coefficient about bbeta(U) in U.fix.1#####
  for(k in 1:P)
  {
    bbeta_true_fix[,k] <- approx(U,DATA[[sn]]$bbeta[,k],xout=U.fix.1,rule=2:1)$y
  }
####extract the estimated Xbeta point to estimate link function
  Xbeta.est=rowSums(bbeta.est[,1:P]*x)
  ######fix the Xbeta point to estimate link function####
  Xbeta.fix<-seq(quantile(Xbeta.est,probs = 0.05),quantile(Xbeta.est,probs = 0.95),diff(c(quantile(Xbeta.est,probs = 0.05),quantile(Xbeta.est,probs = 0.95)))/(Points))
  g.fix.est<-link.kernel(h[2],DATA[[sn]],bbeta.est[,1:P],Para.initial$V,Xbeta.fix)
  link.fix.est <- g.fix.est[,1]
#######estimate the link function under Xbeta.fix####### to evaluate the performance#####
  link_true_fix <- linkFun(Xbeta.fix,case)[[1]]
  print(paste0(mean(abs(g.fix.est[,1]-Xbeta.fix))))  
##########in order to estimate Variance,for binary data, we use p(1-p)########
  link.est<-approx(Xbeta.fix,link.fix.est,xout=DATA[[sn]]$Xbeta,rule=2)$y######## fix link(xbeta) point to estimate Variance
  bbeta.est.P<-matrix(NA,(Points+1),P)
  for(k in 1:P){
#bbeta.est at Points from dp1 in order to plot
    bbeta.est.P[,k]<-approx(U,bbeta.est[,k],xout=dp1,rule=2)$y
    }
  link.est.P<-matrix(NA,Points+1,1)
####link.fix.est at Points from dp2 in order to plot and Variance.fix.est at Points from dp3
  link.est.P<-approx(DATA[[sn]]$Xbeta,link.est,xout=dp2,rule=2:1)$y
  link.fix<-seq(quantile(link.est,probs = 0.2),quantile(link.est,probs = 0.8),diff(c(quantile(link.est,probs = 0.2),quantile(link.est,probs = 0.8)))/Points)
  if(case==4){
    Variance.fix.est <-try(v.hf(h[3],DATA[[sn]],link.est,link.fix)[,1]-link.fix^2,silent = TRUE)
  	if('try-error' %in% class(Variance.fix.est)){
       next
    }else{
      Variance_true_fix <-link_true_fix*(1-link_true_fix) 	
      #Variance.fix.est[Variance.fix.est<0] <- min(Variance.fix.est[Variance.fix.est>0])
      Variance.est <-approx(link.fix.est,Variance.fix.est,xout=DATA[[sn]]$g[,1],rule=2)$y
      Variance.est.P<-approx(DATA[[sn]]$g[,1],Variance.est,xout=dp3,rule=2)$y	
      parameter_est <- cbind(bbeta.fix.1.est,link.fix.est,Variance.fix.est,bbeta_true=bbeta_true_fix,g_true=            link_true_fix,Variance_true =Variance_true_fix)
      parameter_fix <- cbind(bbeta.est.P,link.est.P,Variance.est.P)
    }  
  } else {
    Variance.fix.est <-try(v.hf(h[3],DATA[[sn]],link.est,link.fix)[,1]-link.fix^2,silent = TRUE)
    if('try-error' %in% class(Variance.fix.est)){
       next
    }else{ 
    Variance.fix.est = Variance.fix.est; Variance_true_fix <- link.fix*(1-link.fix)  
    Variance.fix.est[Variance.fix.est<0] <- min(Variance.fix.est[Variance.fix.est>0])
    Variance.est <-approx(link.fix.est,Variance.fix.est,xout=DATA[[sn]]$g[,1],rule=2)$y
    Variance.est.P<-approx(DATA[[sn]]$g[,1],Variance.est,xout=dp3,rule=2)$y	
    len = length(Test_Data$Xbeta)
    g_predict <- approx(DATA[[sn]]$Xbeta,link.est,xout =Test_Data$Xbeta,rule =2)$y
    y_predict <- rep(0,round(N/10))
    for(i in 1:len){y_predict[i] <- ifelse(g_predict[i]>=0.5,1,0)}
      accuracy = Test_Data$Y-y_predict
    predict_accuracy = length(which(accuracy==0))/len
    parameter_est <- cbind(bbeta.fix.1.est,link.fix.est,Variance.fix.est,predict_accuracy =predict_accuracy ,bbeta_true=bbeta_true_fix,link_true= link_true_fix,Variance_true =Variance_true_fix)
      parameter_fix <- cbind(bbeta.est.P,link.est.P,Variance.est.P) 
    }
      
  }   
  return(list(result1 = parameter_est,result.fix = parameter_fix))
}

vary.fit<-function(sn)
{
  colnames(DATA[[sn]]$X) <- c(1:P);U = DATA[[sn]]$U
  nk=nk0
  gg <- gam(Y~s(U,by=X.1,bs="cr",k=nk0+2)+s(U,by=X.2,bs="cr",k=nk0+2),data = data.frame(DATA[[sn]]),family=binomial())
  # predcit type = terms is X*bbeta
  bbeta.est <- predict(gg,type="terms",se.fit=T)$fit/DATA[[sn]]$X;
  #######fix point to estimate
  bbeta.fix.1.est<-matrix(NA,(Points+1),P);bbeta_true_fix <- matrix(NA,(Points+1),P) 
  U.fix.1<-seq(quantile(U,probs = 0.1),quantile(U,probs = 0.9),diff(c(quantile(U,probs = 0.1),quantile(U,probs = 0.9)))/Points)
  for(k in 1:P)
  {
    bbeta.fix.1.est[,k] <- approx(U,bbeta.est[,k],xout=U.fix.1,rule=2)$y
  }
#####Find the true function coefficient about bbeta(U) in U.fix.1#####
  for(k in 1:P)
  {
    bbeta_true_fix[,k] <- approx(U,DATA[[sn]]$bbeta[,k],xout=U.fix.1,rule=2:1)$y
  }
  
  link.est<- predict(gg,type="response",se.fit=T,newdata = data.frame(DATA[[sn]]))$fit 
  ####extract the estimated Xbeta point to estimate link function
  Xbeta.est=rowSums(bbeta.est[,1:P]*DATA[[sn]]$X)
  #############fix the Xbeta point to estimate link function###########
  Xbeta.fix<-seq(quantile(Xbeta.est,probs = 0.05),quantile(Xbeta.est,probs = 0.95),diff(c(quantile(Xbeta.est,probs = 0.05),quantile(Xbeta.est,probs = 0.95)))/Points)
  link_true_fix<-linkFun(Xbeta.fix,case)[[1]]
  link.fix.est<-approx(DATA[[sn]]$Xbeta,link.est,xout=Xbeta.fix,rule=2)$y
  #############fix the link point to estimate variance function##########
  link.fix<-seq(quantile(link.est,probs = 0.05),quantile(link.est,probs = 0.95),diff(c(quantile(link.est,probs = 0.05),quantile(link.est,probs = 0.95)))/Points) 
  Variance.fix.est <- link.fix*(rep(1,Points+1)-link.fix)
  Variance_true_fix <- link_true_fix*(rep(1,Points+1)-link_true_fix)
  Variance_true_fix[Variance_true_fix< 0]<-min(Variance_true_fix[Variance_true_fix>0])
  bbeta.est.P<-matrix(NA,(Points+1),P)
  for(k in 1:P){
#bbeta.est at Points from dp1 in order to plot
      bbeta.est.P[,k]<-approx(U,bbeta.est[,k],xout=dp1,rule=2)$y
  }
  link.est.P<-matrix(NA,Points+1,1)
####link.fix.est at Points from dp2 in order to plot and Variance.fix.est at Points from dp3 
  link.est.P<-approx(Xbeta.fix,link.fix.est,xout=dp2,rule=2)$y 
  Variance.est.P<-approx(link.fix,Variance.fix.est,xout=dp3,rule=2)$y
  if(case==4){
    parameter_est <- cbind(bbeta.fix.1.est,link.fix.est,Variance.fix.est,bbeta_true=bbeta_true_fix,g_true= link_true_fix,Variance_true =Variance_true_fix)
    parameter_fix <- cbind(bbeta.est.P,link.est.P,Variance.est.P)  
  } else {
    len = length(Test_Data$Xbeta)
    g_predict <- approx(DATA[[sn]]$Xbeta,link.est,xout =Test_Data$Xbeta,rule =2)$y
    y_predict <- rep(0,round(N/10))
    for(i in 1:len){y_predict[i] <- ifelse(g_predict[i]>=0.5,1,0)}
      accuracy = Test_Data$Y-y_predict
    predict_accuracy = length(which(accuracy==0))/len
    parameter_est <- cbind(bbeta.fix.1.est,link.fix.est,Variance.fix.est,predict_accuracy =predict_accuracy ,bbeta_true=bbeta_true_fix,link_true= link_true_fix,Variance_true =Variance_true_fix)
    parameter_fix <- cbind(bbeta.est.P,link.est.P,Variance.est.P) 
  }
  return(list(result1 = parameter_est,result.fix = parameter_fix))
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