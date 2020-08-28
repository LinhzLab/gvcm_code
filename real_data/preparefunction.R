library(mgcv)
library(fda)

CV1=function(N,KF){
z=rep(1:KF,ceiling(N/KF))[1:N]
set.seed(44);z=sample(z,N)
mm=list();for (i in 1:KF) mm[[i]]=(1:N)[z==i]
return(mm)}

K<-function(h,u,mu)
{ N=length(u)
  if(hType==1){Knel<-1/h*3/4*(1-((u-mu)/h)^2)* ifelse(abs((u-mu)/h)<=1,1,0)}
  if(hType==2){nh<-quantile(unique(abs(u-mu)[order(abs(u-mu))]),h);Knel<-1/nh*3/4*(1-((u-mu)/nh)^2)* ifelse(abs((u-mu)/nh)<=1,1,0)} 
  return(Knel)
}

#qr(PJ)$rank < min(dim(PJ))
# estimate function

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
v.hf<-function(h,Data,Beta,Beta_dp,g,dp)
{#Data<-DATA[[sn]];Beta<-BETA0;Beta_dp<-BETA0_dp;dp<-dp3;g<-gh[,1]
  y<-Data$Y; x<-Data$X; u<-Data$U; N=length(y); P=ncol(x)
  gf=g;
  yita=(y)^2
  vh=matrix(NA,length(dp),2)
  for(i in 1:length(dp))
  {
    W=c(K(h,g,dp[i])); z<-cbind(rep(1,N),g-dp[i])
    PJ3=t(z)%*%(W*z)#; ei3=eigen(PJ3)$values
    #if(any(!is.finite(PJ3))==TRUE){
    #vh[i,]=c(Beta_dp$V[i],4)}else if('try-error' %in% class(try((sqrt(N)*solve(sqrt(N)*PJ3+diag(10^-5,2))),silent=T))){
    #vh[i,]=c(Beta_dp$V[i],4)}else{
    vh[i,]=(sqrt(N)*solve(sqrt(N)*PJ3+diag(10^-5,2)))%*%t(z)%*%(W*yita)#}
    #ifelse(qr(PJ)$rank==dim(PJ)[1],linkh[i,]<-solve( PJ )%*%t(z)%*%W%*%y ,linkh[i,]<-Beta$g[i,])
  }
  return(vh)
}
calculate_gvcm_parallel<-function(sn) 
{#sn=1
  X<-DATA[[sn]]$X;U<-DATA[[sn]]$U;P=ncol(DATA[[sn]]$X);N=length(DATA[[sn]]$U)
  
  BETA<-list(vc=gm1[[sn]][,1:P],g=cbind(gm1[[sn]][,P+1],gm1[[sn]][,P+2]))
  
  dp1<-seq(range(U)[1],range(U)[2],diff(range(U))/199)
  dp2<-seq(range(rowSums(X*BETA$vc))[1],range(rowSums(X*BETA$vc))[2],diff(range(rowSums(X*BETA$vc)))/199)
  g_dp<-approx(rowSums(X*BETA$vc),BETA$g[,1],xout=dp2)$y;dg_dp<-approx(rowSums(X*BETA$vc),BETA$g[,2],xout=dp2)$y
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
  
  while(tol>10^-4&&num<=200)
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
     vch[,k]<-approx(dp1,vch_dp[,k],xout=U,rule=2)$y;vch[,P+k]<-approx(dp1,vch_dp[,P+k],xout=U,rule=2)$y
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
  Mvh_dp<-v.hf(h[3],DATA[[sn]],BETA1,BETA1_dp,gh[,1],dp3)
  vh0_dp<-Mvh_dp[,1]
  
  #vh_dp<-vh0_dp-dp3^2
  #vh_dp[vh_dp<0]<-min(vh_dp[vh_dp>0])
  #vh<-approx(dp3,vh_dp,xout=gh[,1],rule=2)$y   
  
  vh0<-approx(dp3,vh0_dp,xout=gh[,1],rule=2)$y 
  vh<-vh0-gh[,1]^2
  vh[vh<0]<-min(vh[vh>0])
  #-----------------------------------------------------
  #cat("seed",sn,"del",tol,"iter",num,"\n")
  #par(mfrow=c(2,P))
  #for(p in 1:P){plot(dp1,vch_dp[,p],ylim=c(min(c(vch_dp[,p],DATA[[sn]]$vc[,p])),max(c(vch_dp[,p],DATA[[sn]]$vc[,p]))));lines(U[order(U)],DATA[[sn]]$vc[order(U),p],col="red")};plot(dp2,gh_dp[,1],ylim=c(min(gh_dp[,1],DATA[[sn]]$g[,1]),max(gh_dp[,1],DATA[[sn]]$g[,1])));lines(DATA[[sn]]$SInd[order(DATA[[sn]]$SInd)],DATA[[sn]]$g[order(DATA[[sn]]$SInd),1],col="red")
  #dp2.0=rowSums(X*vch[,1:P])
  #dp2=seq(range(dp2.0)[1],range(dp2.0)[2],diff(range(dp2.0))/199)
  #V=BETA1$g[,1]*(1-BETA1$g[,1]);tor1=10^-4;V[which(V<tor1)]<-tor1;BETA1$V<-V
  
  #gh_dp<-link.hf(h[3],DATA[[sn]],BETA1,BETA1_dp,vch,dp2)
  #gh<-cbind(spline(dp2,gh_dp[,1],xout=dp2.0)$y,spline(dp2,gh_dp[,2],xout=dp2.0)$y)
  
  #BETA1<-list(vc=vch[,1:P],g=gh)
  #beta1<-matrix(c(BETA1$vc,BETA1$g[,1],BETA1$g[,2]),N,P+2)
  #BETA1_dp<-list(vc=vch_dp[,1:P],g=gh_dp)
  #beta1_dp<-matrix(c(BETA1_dp$vc,BETA1_dp$g[,1],BETA1_dp$g[,2]),length(dp1),P+2) 
  BETA1<-list(vc=vch[,1:P],g=gh,v=vh)
  beta1<-matrix(c(BETA1$vc,BETA1$g[,1],BETA1$v),N,P+2)
  
  result1<-matrix(cbind(beta1[,1:(P+2)],rep(tol,N)),nrow=N,ncol=P+3)
  #result1<-matrix(cbind(beta1_dp[,1:(P+1)],rep(tol,nrow(beta1_dp)),dp1,dp2),nrow=nrow(beta1_dp),ncol=P+4)
  
  #res<-list(result=result,result_dp=result_dp)
  return(result1)
}

ini.fit<-function(DATA,nk0,initial)
{#DATA<-d1;nk0=3;P=ncol(DATA$X)
 #id0<-sample(1:N,4000,T,rep(1/N,N));DATA<-d2[id0,];nk0=3;P=ncol(DATA$X)
  x=DATA$X;y=DATA$Y;u<-DATA$U;P=ncol(x);N=length(y)
  nk=nk0
  N=length(y);P=ncol(x);m0<-matrix(c(y,rep(0,N)),N,2);num=1;tol=1
  m<-m0[,1]
  if(P==2){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3),family=binomial())}
  if(P==3){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3),family=binomial())}
  if(P==4){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3),family=binomial())}
  if(P==5){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3)+s(u,by=x[,5],bs="cr",k=nk+3),family=binomial())}
  if(P==6){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3)+s(u,by=x[,5],bs="cr",k=nk+3)+s(u,by=x[,6],bs="cr",k=nk+3),family=binomial())}
  if(P==7){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3)+s(u,by=x[,5],bs="cr",k=nk+3)+s(u,by=x[,6],bs="cr",k=nk+3)+s(u,by=x[,7],bs="cr",k=nk+3),family=binomial())}
  if(P==8){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3)+s(u,by=x[,5],bs="cr",k=nk+3)+s(u,by=x[,6],bs="cr",k=nk+3)+s(u,by=x[,7],bs="cr",k=nk+3)+s(u,by=x[,8],bs="cr",k=nk+3),family=binomial())}
  if(P==9){gg<-gam(m~s(u,by=x[,1],bs="cr",k=nk+3)+s(u,by=x[,2],bs="cr",k=nk+3)+s(u,by=x[,3],bs="cr",k=nk+3)+s(u,by=x[,4],bs="cr",k=nk+3)+s(u,by=x[,5],bs="cr",k=nk+3)+s(u,by=x[,6],bs="cr",k=nk+3)+s(u,by=x[,7],bs="cr",k=nk+3)+s(u,by=x[,8],bs="cr",k=nk+3)+s(u,by=x[,9],bs="cr",k=nk+3),family=binomial())}
  plot(gg,page=1)
  #gg<-dget("~/Documents/gvcm phone/result/initial detect/6C x5large")  
  #dput(gg,"~/Documents/gvcm phone/result/initial detect/6C x5large") 
  g0<-predict(gg,type="iterms",se.fit=T)$fit;pm<-predict(gg,type="link",se.fit=T)$fit
  g0.se<-predict(gg,type="iterms",se.fit=T)$se.fit
  G0=rowSums(g0);bG0<-bsplineS(G0,quantile(G0,seq(0,1,1/nk)),nderiv=0)
  m01=exp(G0)/(1+exp(G0));dm01=exp(G0)/(1+exp(G0))^2
  
  m0<-matrix(c(m01,dm01),N,2);m1<-m0
  
  #par(mfrow=c(2,3))
  for(i in 1:P)
  {
   g0[,i]<-g0[,i]/x[,i];g0[which(x[,i]==0),i]<-g0[which(x[,i]==0),i]
   plot(u,g0[,i])
  }
  plot(G0,m1[,1])
  upband<-g0+1.96*g0.se
  loband<-g0-1.96*g0.se
  
  beta1<-data.frame(g0);beta1$m<-c(m1[,1]);beta1$dm<-m1[,2]
  return(beta1)
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