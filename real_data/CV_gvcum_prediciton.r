rm(list = ls())
setwd("/home/csr2/gvcm_phone")
library(splines)
library(splines2)
library(Matrix)
library(nlme)
library(mgcv)
library(fda)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(modelr)

source('preparefunction_2.R')
#初值估计
ddata<-d9reorder
K_F = 5; nk0 = 5; CV_pre_gm<-list();DATA<-list();gm1<-list()
hType=2;h=c(0.5,0.5,0.5);NORM=3
set.seed(3)
split_K <- crossv_kfold(data.frame(cbind(ddata$X,ddata$U,ddata$Y)),k = K_F);
re = list();Sgm = list(); train_data = list(); test_data = list()
for(k in 1:K_F){ 
 index = split_K$train[[k]]$idx;len <- length(index);index_test = split_K$test[[k]]$idx 
 test_data[[k]] <-list(X=ddata$X[index_test,],U=ddata$U[index_test],Y=ddata$Y[index_test])
 train_data[[k]] <- list(X=ddata$X[index,],U=ddata$U[index],Y=ddata$Y[index])
#interpolate to compute inital parameter
 inital.par <- weight.sp.fit(k,nk0)
 Sgm[[k]] <- data.frame(beta = inital.par[[1]], m =inital.par[[2]], dm = inital.par[[3]])
  print(k)
  k = k+1
 }
DATA = train_data;gm1 <- Sgm;
#dput(gm1,"/home/csr2/gvcm_phone/initial estimater by splines");
cl.cores <- detectCores()
numWorkers=K_F
cl <- makeCluster(numWorkers)
registerDoParallel(cl)
# choose spline as inital parameter into kernel repeat one step.
re <- foreach(iter=1:K_F,.packages=c('foreach','parallel')) %dopar% calculate_gvcm_parallel(iter)
stopCluster(cl)
stopImplicitCluster()
#dput(re,"/home/csr2/gvcm_phone/mcv_5k_seed_3_re_h0.27_norm3_sp")
pred_error_test=matrix(NA,K_F,1);pred_error_train=matrix(NA,K_F,1);P <- ncol(ddata$X)
for(k in 1:K_F){
len1 <- length(test_data[[k]]$U);TBeta <- matrix(NA,len1,P);TG <- matrix(NA,len1,1)
for(p in 1:P){TBeta[,p]<-approx(train_data[[k]]$U,re[[k]][,p],xout=test_data[[k]]$U,rule=2)$y}
TG = approx(rowSums(train_data[[k]]$X*re[[k]][,1:P]),re[[k]][,P+1],xout=rowSums(test_data[[k]]$X*TBeta),rule=2)$y
pred_error_train[k,] = mean(sapply(re[[k]][,P+1]-train_data[[k]]$Y,function(x) abs(x)),na.rm = TRUE)
pred_error_test[k,] = mean(sapply(TG-test_data[[k]]$Y,function(x) abs(x)),na.rm = TRUE)
k = k+1
print(k)
}
avg_pre_error_train= mean(pred_error_train)
avg_pre_error_test = mean(pred_error_test)
avg_pre_error_train
avg_pre_error_test
#save.image(file='gmcv_5k_seed_3.Rdata_h0.27_norm3_sp')

#DATA[[1]] = list(X = ddata$X[1:10000,],U=as.matrix(ddata$U[1:10000]),Y = ddata$Y[1:10000])
#gm0<-ini.fit(DATA[[1]],3,1);dput(gm0,"/home/stu1/gvcm_phone/initial estimater d9reorder");P=ncol(gm0[,1:(ncol(gm0)-2)]);N=nrow(gm0)
#gm0<-dget("/home/csr5/gvcm phone/result/initial estimater d9reorder");P=ncol(gm0[,1:(ncol(gm0)-2)]);N=nrow(gm0)
#par(mfrow=c(2,P+1)); for(i in 1:P){plot(d12$U,gm0[,i])}; plot(rowSums(gm0[,1:P]*d12$X),gm0[,1+P])


