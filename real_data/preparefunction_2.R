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
    if(hType==1){Knel<-1/h*3/4*(1-((u-mu)/h)^2)* ifelse(abs((u-mu)/h)<=1,1,0)}
    if(hType==2){nh<-quantile(unique(abs(u-mu)[order(abs(u-mu))]),h);Knel<-1/nh*3/4*(1-((u-mu)/nh)^2)* ifelse(abs((u-mu)/nh)<=1,1,0)} 
    return(Knel)
  }

  #qr(PJ)$rank < min(dim(PJ))
  
  # estimate \beta function
  vc.hf<-function(h,Data,Beta,Beta_dp,dp)
  { #Data<-DATA[[sn]];Beta<-BETA0;dp<-dp1
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    vc0<-cbind(Beta$vc,matrix(0,N,P)); g<-Beta$g[,1]; dg<-Beta$g[,2]; V<-Beta$V
    vc1=matrix(NA,length(dp),2*P)
    for(i in 1:length(dp))
    {
      #M.U<-diag(u-dp[i]); GMU<-matrix(c(x,M.U%*%x),N,2*P) 
      M.U<-c(u-dp[i]); GMU<-matrix(c(x,M.U*x),N,2*P)
      #W1=diag(K(h,u,dp[i])*dg/V) ; W2=diag(K(h,u,dp[i])*dg^2/V) 
      W1=c(K(h,u,dp[i])*dg/V) ; W2=c(K(h,u,dp[i])*dg^2/V) 
      yita=y-g+dg*c(rowSums(x*vc0[,1:P]))
      #PJ=t(GMU)%*%W2%*%GMU
      PJ=t(GMU)%*%(W2*GMU)
      #print(PJ)
      #if(anyNA(PJ)){
      #vc1[i,]<-Beta_dp$vc[i,]}else if('try-error' %in% class(try(solve(PJ),silent=T))){
      # vc1[i,]<-Beta_dp$vc[i,]}else{
      #vc1[i,]=solve(PJ+diag(10^-5,2*P))%*%t(GMU)%*%W1%*%yita#}  
      vc1[i,]=solve(PJ+diag(10^-5,2*P))%*%t(GMU)%*%(W1*yita)#} 
      #if(NORM==1){Ident=sum((vc1[i,1:P])^2)^(1/2); vc1[i,]=vc1[i,]/Ident} # ||vc(U)||=1 for any U
    }
    # if(NORM==2){Ident=rowSums((vc1[,1:P])^2)[which.min(dp)]^(1/2); vc1=sign(vc1[which.min(dp),1])*vc1/Ident }# ||vc(0)||=1 
    # if(NORM==3){Ident=sum((vc1[,1])^2)^(1/2); vc1=vc1/Ident }# ||vc1(U)||=1 
    #Ident=(sum(colSums(c(diff(range(u))/N)*vc1[,1:P]^2)))^(1/2) ; vc1=vc1/Ident # sum_j(integeral(vc^2))=1
    return(vc1)
  }

  link.hf<-function(h,Data,Beta,Beta_dp,vc,dp)
  {#Data<-DATA[[sn]];Beta<-list(vc=CVgm[[sn]][,1:P],V=CVgm[[sn]][,P+1]*(1-CVgm[[sn]][,P+1]));vc<-Beta$vc;dp<-dp2
   #dp<-seq(range(rowSums(Data$X*Beta$vc))[1],range(rowSums(Data$X*Beta$vc))[2],diff(range(rowSums(Data$X*Beta$vc)))/199)
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    vcf=vc[,1:P]; V<-Beta$V
    linkh=matrix(NA,length(dp),2)
   for(i in 1:length(dp))
    {
     # W=diag(K(h,rowSums(x*vcf),dp[i])/V); z<-cbind(rep(1,N),rowSums(x*vcf)-dp[i])
     # PJ=t(z)%*%W%*%z
      W=c(K(h,rowSums(x*vcf),dp[i])/V); z<-cbind(rep(1,N),rowSums(x*vcf)-dp[i])
      PJ=t(z)%*%(W*z)
      # if(qr(PJ)$rank<min(dim(PJ))){
      #  linkh[i,]<-Beta_dp$g[i,]}else if('try-error' %in% class(try(solve(PJ),silent=T))){
      #   linkh[i,]<-Beta_dp$g[i,]}else{
      #linkh[i,]=solve(PJ+diag(10^-5,2))%*%t(z)%*%W%*%y#}  
      linkh[i,]=solve(PJ+diag(10^-5,2))%*%t(z)%*%(W*y)#}
      #ifelse(qr(PJ)$rank==dim(PJ)[1],linkh[i,]<-solve( PJ )%*%t(z)%*%W%*%y ,linkh[i,]<-Beta$g[i,])
    }
    return(linkh)
  }
  
v.hf<-function(h,Data,g,dp)
  {#Data<-DATA[[sn]];Beta<-BETA0;Beta_dp<-BETA0_dp;dp<-dp3;g<-gh[,1]
    y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
    yita=(y)^2
    vh=matrix(NA,length(dp),2)
    for(i in 1:length(dp))
    {
      W=c(K(h,g,dp[i])); z<-cbind(rep(1,N),g-dp[i])
      PJ3=t(z)%*%(W*z); vh[i,]=(sqrt(N)*solve(sqrt(N)*PJ3+diag(10^-5,2)))%*%t(z)%*%(W*yita)
    }
    return(vh)
  }

#use spline method to get the inital parameters of \beta, \dot g and  g; P: dimension of X
  spline.inital<-function(Bu,y,x,P,sn,nk)
  {
#objective function, nk: the number of basis function, thus total parameter (p+1)*nk   	
    objec_f <- function(parameter){
      alpha_2 <- matrix(parameter[-(1:nk)],nrow = P,byrow = T)
# get x\beta(u)      
      xu <- rowSums(x*(Bu%*%t(alpha_2)))
      Bxu_1 <- bSpline(xu,df = nk,intercept = FALSE)
      expression_1 <- y-Bxu_1%*%parameter[1:nk]
      expression <- colSums(expression_1*expression_1)
      }
    inital_parameter <- rep(0.1,((P+1)*nk))
    parameter <- optim(inital_parameter,objec_f)$par
    beta_hat <- Bu%*%t(matrix(parameter[-(1:nk)],nrow = P,byrow = T))
    xu <- rowSums(x*beta_hat); Bxu_1 <- bs(xu,df = nk,intercept = FALSE)
    g_hat <-  Bxu_1%*%parameter[1:nk]
    return(list(beta_hat,g_hat,parameter))
  }
  #spline.inital<- spline.inital(sn,3)
  #mu <- spline.inital[[2]]
  sp_Variance <- function(y,sn,nk,mu){
    B_base <- bSpline(mu,df = nk,intercept = FALSE)
    objec_f <- function(parameter){
    	expression_0 <- y*y-B_base%*%parameter
      expression <- colSums(expression_0*expression_0)
      }
   in_parameter <- rep(0.1,nk)
   parameter <- optim(in_parameter,objec_f)$par;variance <- B_base%*%parameter 
   return(variance)
  }
#derive g and \dot g function by (2.6) in paper 
  link.kernel<-function(h,beta_hat,dp,Vari,y,x,N){
    linkh=matrix(NA,length(dp),2)
    for(i in 1:length(dp)){
      W=c(K(h,rowSums(x*beta_hat),dp[i])/(Vari+rep(0.00001,N)));z<-cbind(rep(1,N),rowSums(x*beta_hat)-dp[i])
      PJ=t(z)%*%(z*W)
      linkh[i,]=solve(PJ+diag(10^-5,2))%*%t(z)%*%(W*y)#}
      #print(linkh[i,])
      i = i+1
      }
      return(linkh)
  }

#Vari <- sp_Variance(sn,nk0,mu)
# estimate inital parameter under weight QMLE with spline method  
 weight.sp.fit <- function(sn,nk0){
    Data = train_data[[sn]]; y<-Data$Y; x <- Data$X; u<- Data$U; N=length(y); P=ncol(x)	
    nk = nk0; Bu <- bSpline(u,df = nk,intercept = FALSE)
    sp.inital<- spline.inital(Bu,y,x,P,sn,nk)
    mu <- sp.inital[[2]];Vari <- sp_Variance(y,sn,nk,mu);Vari <- Vari-mu^2
    op.f <- function(parameter){
      alpha_2 <- matrix(parameter[-(1:nk)],nrow = P,byrow = T)
      xu <- rowSums(x*(Bu%*%t(alpha_2)))
      Bxu_1 <- bSpline(xu,df = nk,intercept = FALSE)
      expression_1 <- y-Bxu_1%*%parameter[1:nk]
      expression <- colSums((expression_1*expression_1)/(Vari+rep(0.00001,N)))
      return(expression)
      }
    parameter <- optim(sp.inital[[3]],op.f)$par
    beta_hat <- Bu%*%t(matrix(parameter[-(1:nk)],nrow = P,byrow = T))
    xu <- rowSums(x*beta_hat); Bxu_1 <- bSpline(xu,df = nk,intercept = FALSE)
    g_hat <- Bxu_1%*%parameter[1:nk]
    dp1.0=rowSums(x*beta_hat[,1:P])
    dp=seq(range(dp1.0)[1],range(dp1.0)[2],diff(range(dp1.0))/199)
    gdp.deriv = link.kernel(h[2],beta_hat,dp,Vari,y,x,N)[,2]
    g.deriv = approx(dp,gdp.deriv,rule = 1:2,xout=xu)$y
    return(list(beta_hat,g_hat,g.deriv))
    } 
  calculate_gvcm_parallel<-function(sn) 
  {
    X<-DATA[[sn]]$X;U<-DATA[[sn]]$U;P=ncol(DATA[[sn]]$X);N=length(DATA[[sn]]$U)
    
    BETA<-list(vc=gm1[[sn]][,1:P],g=cbind(gm1[[sn]][,P+1],gm1[[sn]][,P+2]))
    dp1<-seq(range(U)[1],range(U)[2],diff(range(U))/199)
    dp2<-seq(range(rowSums(X*BETA$vc))[1],range(rowSums(X*BETA$vc))[2],diff(range(rowSums(X*BETA$vc)))/199)
#choose interpolate points for g and deriv-g   
    g_dp<-approx(rowSums(X*BETA$vc),BETA$g[,1],rule = 1:2,xout=dp2)$y
    dg_dp<-approx(rowSums(X*BETA$vc),BETA$g[,2],rule = 1:2,xout=dp2)$y
    link_dp<-cbind(g_dp,dg_dp)    
    vc_dp<-matrix(NA,200,P);dvc_dp<-matrix(NA,200,P)
    for(k in 1:P)
    {
     vc_dp[,k]<-approx(U,BETA$vc[,k],xout=dp1)$y#;dvc_dp[,k]<-approx(U,BETA$dvc[,k],xout=dp1)$y   
    }
    BETA_dp<-list(vc=vc_dp,g=link_dp)
    BETA1<-BETA
    BETA1_dp<-BETA_dp
    num=1;tol=1
    while(tol>10^-2&&num<=1)
    {
      BETA0<-BETA1
      BETA0_dp<-BETA1_dp
      V=BETA0$g[,1]*(1-BETA0$g[,1]);tor1=10^-4;V[which(V<tor1)]<-tor1;BETA0$V<-V
      beta0_dp<-matrix(cbind(BETA0_dp$vc,BETA0_dp$g[,1],BETA0_dp$g[,2]),length(dp1),P+2)
      beta0<-matrix(cbind(BETA0$vc,BETA0$g[,1],BETA0$g[,2]),N,P+2)      
      vch_dp<-vc.hf(h[1],DATA[[sn]],BETA0,BETA0_dp,dp1)#;  vch_dpt<-vch_dp;  par(mfrow=c(2,P));for(p in 1:P){plot(dp1,vch_dp0[,p])}
      if(NORM==1){for(i in 1:200){Ident=sum((vch_dp[i,1:P])^2)^(1/2); vch_dp[i,]=vch_dp[i,]/Ident}} # ||vc(U)||=1 for any U
      if(NORM==1.5){for(i in 1:200){Ident=sum((vch_dp[i,1:P])^2)^(1/2); vch_dp[i,]=sign(vch_dp[i,1])*vch_dp[i,]/Ident}} # ||vc(U)||=1 for any U
      if(NORM==4){for(i in 1:200){Ident=sum((vch_dp[i,1:P])^2)^(1/2); vch_dp[i,]=sign(vch_dp[i,1])*vch_dp[i,]/Ident}} # ||vc(U)||=1 for any U
      if(NORM==2){Ident=rowSums((vch_dp[,1:P])^2)[which.min(U)]^(1/2); vch_dp=sign(vch_dp[which.min(U),1])*vch_dp/Ident} # ||vc(0)||=1 
      if(NORM==3){Ident=vch_dp[which(abs(U-median(U))==min(abs(U-median(U))))[1],1]; vch_dp=vch_dp/Ident} # vc1(0)=1 
      vch<-matrix(NA,N,2*P)
      for(k in 1:P)
      {
       vch[,k]<-approx(dp1,vch_dp[,k],xout=U,rule=1:2)$y;vch[,P+k]<-approx(dp1,vch_dp[,P+k],xout=U,rule=2)$y
      }    
     
      dp2.0=rowSums(X*vch[,1:P])
      dp2=seq(range(dp2.0)[1],range(dp2.0)[2],diff(range(dp2.0))/199)
      
      gh_dp<-link.hf(h[2],DATA[[sn]],BETA0,BETA0_dp,vch,dp2)
      gh<-cbind(spline(dp2,gh_dp[,1],xout=dp2.0)$y,spline(dp2,gh_dp[,2],xout=dp2.0)$y)
      
      BETA1<-list(vc=vch[,1:P],g=gh)
      beta1<-matrix(c(BETA1$vc,BETA1$g[,1],BETA1$g[,2]),N,P+2)
      BETA1_dp<-list(vc=vch_dp[,1:P],g=gh_dp)
      beta1_dp<-matrix(c(BETA1_dp$vc,BETA1_dp$g[,1],BETA1_dp$g[,2]),length(dp1),P+2)
      
      tol=max(abs(beta0_dp-beta1_dp)[,1:(P+1)])
      print(paste("seed",sn,"del",tol,"iter",num,seq=""))
      #par(mfrow=c(2,P));for(p in 1:P){plot(dp1,vch_dp[,p])};plot(dp2,gh_dp[,1])
      num = num + 1
    }#end while
#-----------variance function------------------------
    dp3=seq(range(gh[,1])[1],range(gh[,1])[2],diff(range(gh[,1]))/199)
    Mvh_dp<-v.hf(h[3],DATA[[sn]],gh[,1],dp3)
#derive the variance at dp3    
    vh0_dp<-Mvh_dp[,1]
    
    #vh_dp<-vh0_dp-dp3^2
    #vh_dp[vh_dp<0]<-min(vh_dp[vh_dp>0])
    #vh<-approx(dp3,vh_dp,xout=gh[,1],rule=2)$y   
    
    vh0<-approx(dp3,vh0_dp,xout=gh[,1],rule=2)$y 
    vh<-vh0-gh[,1]^2
    vh[vh<0]<-min(vh[vh>0])
    BETA1<-list(vc=vch[,1:P],g=gh,v=vh)
    beta1<-matrix(c(BETA1$vc,BETA1$g[,1],BETA1$v),N,P+2)
    result1<-matrix(cbind(beta1[,1:(P+2)],rep(tol,N)),nrow=N,ncol=P+3)
    return(result1)
  }

  ##############################
  d1<-dget("phoneData");do<-dget("phoneDatao");
  #d2=list(X=scale(d1$X[,-6]),U=d1$U,Y=d1$Y)
  #d3=list(X=scale(d1$X[,-6]),U=c(scale(d1$U)),Y=d1$Y)
  #d4=list(X=scale(d1$X),U=d1$U,Y=d1$Y)
  #d5=list(X=scale(d1$X),U=scale(d1$U),Y=d1$Y)
  #d6=list(X=scale(d1$X[,-2]),U=scale(d1$U),Y=d1$Y)
  #d7=list(X=scale(d1$X[,-5]),U=scale(d1$U),Y=d1$Y)
  #d8=list(X=scale(cbind(d1$X[,-1],d1$U)),U=scale(d1$X[,1]),Y=d1$Y)
  d9=list(X=scale(cbind(d1$X[,-c(1,2)],d1$U)),U=scale(d1$X[,1]),Y=d1$Y)
  d9reorder=list(X=cbind(d9$X[,2],d9$X[,1],d9$X[,3:5]),U=d9$U,Y=d9$Y)
  d09=list(X=scale(cbind(d1$X[,-c(1,2)],d1$U)),U=scale(log(d1$X[,1])),Y=d1$Y)
  d1$U<-d1$U/1000;d1$U[which(d1$U==0)]<-1;
  d009=list(X=scale(cbind(d1$X[,-c(1,2)],log(d1$U))),U=scale(log(d1$X[,1]/1000)),Y=d1$Y);
  #d10=list(X=scale(cbind(d1$X[,-c(3)],d1$U)),U=scale(d1$X[,3]),Y=d1$Y)
  #d11=list(X=scale(cbind(d1$X[,-c(2,3)],d1$U)),U=scale(d1$X[,3]),Y=d1$Y)
  #d12=list(X=scale(cbind(d1$X[,-c(2,3,4)],d1$U)),U=scale(d1$X[,3]),Y=d1$Y)