rm(list = ls())
#setwd("/home/csr2/gvcmsim")
#setwd("/data/liujx/gvcmsim")
setwd("/home/LiuJX/gvcmsim")
library(snow)
library(snowfall)
library(splines2)
library(splines)
library(Matrix)
library(fda)
library(MASS)
source("pre_sim_3_case8.R")
source("generate_data_3.R")

Points = 199;N=8000; NORM = 3;S= 500;nk2 = 4;P = 2;Sigma<-diag(P);nk1 = 5;nk =4;h=c(0.25,0.75,0.45)
dim1<-c(paste("seed",1:S));dim3<-c("beta1","beta2","g",'Variance',"beta1_true","beta2_true","g_true",'Variance_true')
dim3.fix <- c("beta.1.fix","beta.2.fix","g.fix",'Variance.fix')
result1<-array(NA,dim=c(S,(Points+1),length(dim3)),dimnames=list(dim1,NULL,dim3)) 
result.fix<-array(NA,dim=c(S,(Points+1),length(dim3.fix)),dimnames=list(dim1,NULL,dim3.fix));
List.Kern<-list(N.8000=NA,N.15000=NA,N.20000=NA)
sfInit(parallel = TRUE, cpus =25)
while(N<20001){
#####fix point for U:dp1, betaX:dp2, g(betaX):dp3 to estimate#####
set.seed(1)
U<-runif(N,0,1)
dp1<-seq(0,1,diff(c(0,1))/Points)
x<-mvrnorm(N, rep(0, P), Sigma);bbeta=matrix(cbind(sin(0.5*pi*U),cos(0.5*pi*U)),N,P);Xbeta=rowSums(x*bbeta)
dp2<-seq(quantile(Xbeta,probs = 0.05),quantile(Xbeta,probs = 0.95),diff(c(quantile(Xbeta,probs = 0.05),quantile(Xbeta,probs = 0.95)))/Points); g<-linkFun(Xbeta)[[1]]
dp3<-seq(quantile(g,probs = 0.05),quantile(g,probs = 0.95),diff(c(quantile(g,probs = 0.05),quantile(g,probs = 0.95)))/Points)
  
######################################################	
set.seed(210)
test_covariate <-genInput(1000,P);Test_Data = genData(test_covariate)
DATA<-list();input_covariate<-list()
for(i in 1:S){ 
  set.seed(i)
  input_covariate[[i]]<-genInput(N,P)
  DATA[[i]]<-genData(input_covariate[[i]])
}
#cl.cores <- detectCores()
sfLibrary(splines2)
sfLibrary(splines)     # 载入依赖R包MASS
sfExport("DATA",'nk','nk1','nk2','Points',"P","N",'h','NORM')
sfSource('pre_sim_3_case8.R')
initial.estimate <- sfLapply(1:S, spline.initial) 
#path <- paste0("/home/csr2/gvcmsim/initial_estimate_example_3",N,".Rdata")
#save(initial.estimate,file=path)
#load(path) 
gm1 = list();gm1 <- initial.estimate 
sfExport("gm1","DATA","Test_Data",'nk','nk1',"nk2","P","N",'h','Points','dp1','dp2','dp3')
sfSource('pre_sim_3_case8.R')
sfSource('generate_data_3.R')
set.seed(1)
res <- sfLapply(1:S, calculate_gvcm_parallel) 
  #res <- foreach(iter=1:S) %do% calculate_gvcm_parallel(iter)
result1<-array(NA,dim=c(S,(Points+1),length(dim3)),dimnames=list(dim1,NULL,dim3)) 
result.fix<-array(NA,dim=c(S,(Points+1),length(dim3.fix)),dimnames=list(dim1,NULL,dim3.fix)) 
Y_pred<-array(NA,dim=c(S,1),dimnames=list(dim1,"MSE")) 
for(i in 1:S){result1[i,,]<-res[[i]]$result1;result.fix[i,,]<-res[[i]]$result.fix;Y_pred[i]<-res[[i]]$MSE;i=i+1}
if(N==20000){List.Kern$N.20000<-list(result1 = result1,result.fix =result.fix,Y_pred);N=N+10;print(N)}
if(N==15000){List.Kern$N.15000<-list(result1 = result1,result.fix =result.fix,Y_pred);N=N+5000;nk1 = 5;nk =4;h=c(0.25,0.5,0.30);print(N)}
if(N==8000){List.Kern$N.8000<-list(result1 = result1,result.fix =result.fix,Y_pred);N=N+7000;nk1 =5;nk =4;h=c(0.25,0.5,0.38);print(N)}
#if(N==100){List.Kern$N.100<-list(result1 = result1,result.fix =result.fix);N=N+100;print(N)} 
}
sfStop()
#path <- paste0("/home/csr2/gvcmsim/gvcm_example3_300_8000.Rdata")
#path <- paste0("/home/liujx/gvcmsim/gvcm_example1_case_2.Rdata")
path <- paste0("/home/LiuJX/gvcmsim/gvcm_example3_case8.Rdata")
save.image(file=path)
#load(path)  
#source("MSE for One step.R")
#MSE
#MSE1


