rm(list = ls())
#setwd("/home/stu1/gvcm_phone")
setwd("/home/csr2/gvcm_phone")
library(splines)
library(Matrix)
library(nlme)
library(mgcv)
library(fda)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(modelr)
source('preparefunction.R')
#初值估计
ddata<-d9reorder
P=ncol(ddata$X)
K_F = 5; CV_pre_gm<-list();DATA<-list();gm1<-list()
hType=2;h=c(0.5,0.5,0.5);NORM=2
set.seed(5)
split_K <- crossv_kfold(data.frame(cbind(ddata$X,ddata$U,ddata$Y)),k = K_F);
re = list();Sgm = list(); train_data = list(); test_data = list()
for(k in 1:K_F){ 
 index = split_K$train[[k]]$idx;len <- length(index);index_test = split_K$test[[k]]$idx 
 test_data[[k]] <-list(X=ddata$X[index_test,],U=ddata$U[index_test],Y=ddata$Y[index_test])
 colnames(test_data[[k]]$X) = c(1:P)
 train_data[[k]] <- list(X=ddata$X[index,],U=ddata$U[index],Y=ddata$Y[index])
 colnames(train_data[[k]]$X) = c(1:P)
}
vary.fit<-function(sn,nk0)
{
  colnames(DATA[[sn]]$X) <- c(1:P)
  nk=nk0
  if(P==5){gg<-gam(Y~s(U,by=X.1,bs="cr",k=nk+2)+s(U,by=X.2,bs="cr",k=nk+2),data = data.frame(train_data[[sn]]),family=binomial())}
  g_f<- predict(gg,type="link",se.fit=T)$fit;m0 <- sapply(g_f,P_prob);
  pred <- predict(gg,type="response",se.fit=T,newdata = data.frame(test_data[[sn]]))$fit
  return(list(m0,pred))
}
#str(vary.fit(1,3))
nk0 = 3;re_lg = list()
cl.cores <- detectCores()
numWorkers=K_F
cl <- makeCluster(numWorkers)
registerDoParallel(cl)
re_lg <- foreach(iter=1:K_F,.export = c('train_data'),.packages=c('foreach','parallel','mgcv'))%dopar% vary.fit(iter,nk0)
stopCluster(cl)
stopImplicitCluster()
dput(re_lg,"/home/csr2/gvcm_phone/logistic_k5_seed_5")
pred_error_test=matrix(NA,K_F,1);pred_error_train=matrix(NA,K_F,1);
for(k in 1:K_F){
#len1 <- length(test_data[[k]]$U);TBeta <- matrix(NA,len1,P);TG <- matrix(NA,len1,1)
#for(p in 1:P){TBeta[,p]<-approx(train_data[[k]]$U,re_lg[[k]][,p],xout=test_data[[k]]$U,rule=2)$y}
#TG = approx(rowSums(re_lg[[k]][,1:P]*train_data[[k]]$X),re_lg[[k]][,P+1],xout=rowSums(test_data[[k]]$X*test_data[[k]]$U),rule=2)$y	
re = re_lg[[k]]
pred_error_train[k,] = mean(sapply(re[[1]]-train_data[[k]]$Y,function(x) abs(x)),na.rm = T)
pred_error_test[k,] = mean(sapply(re[[2]]-test_data[[k]]$Y,function(x) abs(x)),na.rm = T)
k = k+1
}
avg_pre_error_train= mean(pred_error_train)
avg_pre_error_train
avg_pre_error_test = mean(pred_error_test)
avg_pre_error_test
save.image(file='logistic_k5_seed_5.Rdata')
#re = dget("/home/csr2/gvcm_phone/CV_re_lg_1")
