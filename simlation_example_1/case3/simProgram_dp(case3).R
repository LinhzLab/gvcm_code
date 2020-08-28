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
source("generate_data_1.R")
source("pre_sim_1.R")

Points = 99;case = 3;N=100; NORM = 3;S= 1000;nk=4;nk1= 6;h=c(0.1,0.25,0.7);P = 3;Sigma<-diag(P)
dim1<-c(paste("seed",1:S));dim3<-c("beta1","beta2","beta3","V","beta1_true","beta2_true","beta3_true",'Variance_true')
dim3.fix <- c("beta.1.fix","beta.2.fix","beta.3.fix","g.fix","Variance.fix")
dim2 <- c("g","g_true")
List.Kern<-list(N.100=NA,N.200=NA,N.400=NA)
sfInit(parallel = TRUE, cpus = 10)
while(N<401){
#####fix point for U:dp1, betaX:dp2, g(betaX):dp3 to estimate#####
set.seed(1)
U<-runif(N,0,1)
dp1<-seq(0,1,1/Points)
x<-mvrnorm(N, rep(0, P), Sigma);bbeta=matrix(cbind(U^2+1,cos(pi*U)*cos(pi*U)+rep(0.5,N),2*sin(pi*U)*sin(pi*U)+rep(0.5,N)),N,P);Xbeta=rowSums(x*bbeta)
dp2<-seq(quantile(Xbeta,probs = 0.01),quantile(Xbeta,probs = 0.99),diff(c(quantile(Xbeta,probs = 0.01),quantile(Xbeta,probs = 0.99)))/Points); g<-linkFun(Xbeta,case)[[1]]
dp3<-seq(-1,1,2/Points)

######################################################	
DATA<-list();input_covariate<-list();gm1 = list()
for(i in 1:S){ 
  set.seed(i)
  input_covariate[[i]]<-genInput(N,P)
  DATA[[i]]<-genData(input_covariate[[i]])
}
#cl.cores <- detectCores()
sfLibrary(splines2)
sfLibrary(splines)     # 载入依赖R包MASS
sfExport("DATA",'nk','nk1','Points',"P","N",'h','NORM')
sfSource('pre_sim_1.R')
initial.estimate <- sfLapply(1:S, spline.initial) 
#path <- paste0("/home/csr2/gvcmsim/initial_estimate_example_1",N,".Rdata")
#save(initial.estimate,file=path)
#load(path) 
gm1 <- initial.estimate 
sfExport("case","gm1","DATA",'nk','nk1',"P","N",'h','Points','dp1','dp2','dp3')
sfSource('pre_sim_1.R')
sfSource('generate_data_1.R')
set.seed(1)
res <- sfLapply(1:S, calculate_gvcm_parallel) 
  #res <- foreach(iter=1:S) %do% calculate_gvcm_parallel(iter)
result_g<-array(NA,dim=c(S,(Points+101),length(dim2)),dimnames=list(dim1,NULL,dim2)) 
result_bbeta_var<-array(NA,dim=c(S,(Points+1),length(dim3)),dimnames=list(dim1,NULL,dim3)) 
result.fix<-array(NA,dim=c(S,(Points+1),length(dim3.fix)),dimnames=list(dim1,NULL,dim3.fix)) 
for(i in 1:S){result_g[i,,]<-res[[i]]$result_g;result_bbeta_var[i,,] <- res[[i]]$result_bbeta_var;result.fix[i,,]<-res[[i]]$result.fix;i=i+1} 
if(N==4000){List.Kern$N.4000<-list(result_bbeta_var = result_bbeta_var,result_g = result_g,result.fix =result.fix);N=N+5;print(N)}
if(N==400){List.Kern$N.400<-list(result_bbeta_var = result_bbeta_var,result_g = result_g,result.fix =result.fix);N=N+3600;print(N);nk=5;nk1= 5;h=c(0.05,0.25,1)}
if(N==200){List.Kern$N.200<-list(result_bbeta_var = result_bbeta_var,result_g = result_g,result.fix =result.fix);N=N+200;print(N);nk=5;nk1= 6;h=c(0.1,0.25,0.7)}
if(N==100){List.Kern$N.100<-list(result_bbeta_var = result_bbeta_var,result_g = result_g,result.fix =result.fix);N=N+100;print(N);nk=4;nk1= 6;h=c(0.1,0.25,0.7)} 
}
sfStop()
#path <- paste0("/home/csr2/gvcmsim/gvcm_example1_case_3.Rdata")
#path <- paste0("/home/liujx/gvcmsim/gvcm_example1_case_3.Rdata")
path <- paste0("/home/LiuJX/gvcmsim/gvcm_example1_case_3.Rdata")
save.image(file=path)
#load(path)  
#source("MSE for One step.R")
#MSE
#MSE1


