rm(list = ls())
#setwd("/home/csr2/gvcmsim")
#setwd("/home/liujx/gvcmsim")
setwd("/home/LiuJX/gvcmsim")
library(snow)
library(snowfall)
library(splines2)
library(splines)
library(Matrix)
library(fda)
library(MASS)
source("generate_data_2.R")

Points = 199;N=800; NORM = 3;case= 4; S= 500;nk=6;nk1=5;nk2 = 6;P = 2;Sigma<-diag(P);h=c(0.49,1.8,0.1,0.5);
dim1<-c(paste("seed",1:S))
if(case ==4){
  dim3<-c("beta1","beta2","g","Variance","beta1_true","beta2_true","g_true",'Variance_true')
  List.Kern<-list(N.800=NA,N.1100=NA,N.1500=NA,N.2000=NA);Number =2001
} else { dim3<-c("beta1","beta2","g","Variance","predict_accuracy","beta1_true","beta2_true","g_true",'Variance_true')
  List.Kern<-list(N.800=NA,N.3000=NA,N.10000=NA);Number = 10001
} 
dim3.fix <- c("beta.1.fix","beta.2.fix","g.fix",'Variance.fix')
result1<-array(NA,dim=c(S,(Points+1),length(dim3)),dimnames=list(dim1,NULL,dim3)) 
result.fix<-array(NA,dim=c(S,(Points+1),length(dim3.fix)),dimnames=list(dim1,NULL,dim3.fix));
sfInit(parallel = TRUE, cpus =20)
while(N<Number){
#####fix point for U:dp1, betaX:dp2, g(betaX):dp3 to estimate#####
set.seed(1)
U<-runif(N,0,1)
dp1<-seq(0,1,diff(c(0,1))/Points)
x<-mvrnorm(N, rep(0, P), Sigma);bbeta=matrix(cbind(sin(pi*U),cos(pi*U)),N,P);Xbeta=rowSums(x*bbeta)
dp2<-seq(quantile(Xbeta,probs = 0.05),quantile(Xbeta,probs = 0.95),diff(c(quantile(Xbeta,probs = 0.05),quantile(Xbeta,probs = 0.95)))/Points); g<-linkFun(Xbeta,case)[[1]]
dp3<-seq(quantile(g,probs = 0.15),quantile(g,probs = 0.85),diff(c(quantile(g,probs = 0.15),quantile(g,probs = 0.85)))/Points)

######################################################  
DATA<-list();input_covariate<-list();set.seed(600)
if(case ==5){
  test_covariate <-genInput(N/10,P);Test_Data = genData(test_covariate,case)
} else {Test_Data=list()} 
for(i in 1:S){ 
  set.seed(i)
  input_covariate[[i]]<-genInput(N,P)
  DATA[[i]]<-genData(input_covariate[[i]],case)
}

#cl.cores <- detectCores()
sfLibrary(splines2)
sfLibrary(splines)     # 载入依赖R包MASS
sfExport("DATA",'nk','nk1','nk2','Points',"P","N",'h','NORM','case')
sfSource('pre_sim_2_dh.R')
set.seed(1)
initial.estimate <- sfLapply(1:S, spline.initial) 
#path <- paste0("/home/csr2/gvcmsim/initial_estimate_example_3",N,".Rdata")
#save(initial.estimate,file=path)
#load(path) 
gm1 = list();gm1 <- initial.estimate
sfSource('pre_sim_2_dh.R')
sfSource('generate_data_2.R') 
if(case ==4){
sfExport("gm1","DATA",'nk','nk1','nk2',"P","N",'case','h','Points','dp1','dp2','dp3')
set.seed(1)
res <- sfLapply(1:S, calculate_gvcm_parallel) 
} else { 
sfExport("gm1","DATA","P","N",'case','h','Test_Data','Points','dp1','dp2','dp3')
set.seed(1)
res <- sfLapply(1:S, calculate_gvcm_parallel) 
}
  #res <- foreach(iter=1:S) %do% calculate_gvcm_parallel(iter)
result1<-array(NA,dim=c(S,(Points+1),length(dim3)),dimnames=list(dim1,NULL,dim3)) 
result.fix<-array(NA,dim=c(S,(Points+1),length(dim3.fix)),dimnames=list(dim1,NULL,dim3.fix)) 
for(i in 1:S){result1[i,,]<-res[[i]]$result1;result.fix[i,,]<-res[[i]]$result.fix;i=i+1}
if(case==4){
  if(N==2000){List.Kern$N.2000<-list(result1 = result1,result.fix =result.fix);N=N+10;print(N)}
  if(N==1500){List.Kern$N.1500<-list(result1 = result1,result.fix =result.fix);N=N+500;nk=5;nk1=4;h=c(0.43,1.98,0.06,0.49);print(N)}
  if(N==1100){List.Kern$N.1100<-list(result1 = result1,result.fix =result.fix);N=N+400;nk=5;nk1=4;h=c(0.48,1.98,0.1,0.5);print(N)}
  if(N==800){List.Kern$N.800<-list(result1 = result1,result.fix =result.fix);h=c(0.48,1.8,0.1,0.5);N=N+300;print(N)}
} else {
  if(N==10000){List.Kern$N.10000<-list(result1 = result1,result.fix =result.fix);N=N+10;print(N)}
  if(N==3000){List.Kern$N.3000<-list(result1 = result1,result.fix =result.fix);N=N+7000;h=c(0.35,3,1);print(N)}
  if(N==800){List.Kern$N.800<-list(result1 = result1,result.fix =result.fix);h=c(0.3,2,0.5);N=N+2200;print(N)}
}  
#if(N==100){List.Kern$N.100<-list(result1 = result1,result.fix =result.fix);N=N+100;print(N)} 
}
sfStop()
#path <- paste0("/home/LiuJX/gvcmsim/ex2gvcm_vs_zhang_dh.Rdata")
#save.image(file=path)
#load(path)  
#source("MSE for One step.R")
#MSE
#MSE1


