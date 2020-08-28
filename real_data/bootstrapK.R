rm(list = ls())
setwd("/home/csr5/gvcm phone")
library(mgcv)
library(fda)
library(foreach)
library(doParallel)
library(iterators)

source('preparefunction.R')
#source('preparefunctionS1.R')

#初值估计
ddata<-d9reorder
gm0<-ini.fit(ddata,3,1);dput(gm0,"/home/csr5/gvcm phone/result/initial estimater d9reorder");P=ncol(gm0[,1:(ncol(gm0)-2)]);N=nrow(gm0)
#模型估计
hType=2;h=c(0.5,0.5,0.5);NORM=1
DATA<-list();gm1<-list()
DATA[[1]]<-ddata;gm1[[1]]<-data.frame(gm0)
K.result.scaled<-calculate_gvcm_parallel(1)
#par(mfrow=c(2,P+1)); for(i in 1:P){plot(DATA[[1]]$U,K.result.scaled[,i])}; plot(rowSums(K.result.scaled[,1:P]*DATA[[1]]$X),K.result.scaled[,1+P])
dput(K.result.scaled,"/home/csr5/gvcm phone/result/K estimater d9reorder norm1 55"
######################### Boot-strap 200###################################

S=200;nDP=200
BSd<-ddata;N=nrow(BSd$U);P=ncol(gm0[,1:(ncol(gm0)-2)])

boos<-array(NA,dim=c(N,P+3,S));data3<-array(NA,dim=c(N,P+2,S));bsdata<-array(NA,dim=c(N,P+2,S))
sde<-matrix(NA,N,P+3)

BSdata<-list();BSgm<-list()
for(k in 1:S)
{ 
  set.seed(k)
  index<-sample(1:N,N,T,rep(1/N,N))
  BSdata[[k]]<-list(X=BSd$X[index,],U=BSd$U[index],Y=BSd$Y[index])
  BSvch<-matrix(NA,length(BSdata[[k]]$U),P);BSgh<-matrix(NA,length(BSdata[[k]]$U),2)
  for(p in 1:P){BSvch[,p]<-approx(BSd$U,gm0[,p],xout=BSdata[[k]]$U,rule=2)$y}
  BSgh[,1]<-approx(rowSums(BSd$X*gm0[,1:P]),gm0[,P+1],xout=rowSums(BSdata[[k]]$X*BSvch),rule=2)$y
  BSgh[,2]<-approx(rowSums(BSd$X*gm0[,1:P]),gm0[,P+2],xout=rowSums(BSdata[[k]]$X*BSvch),rule=2)$y
  BSgm[[k]]<-cbind(BSvch,BSgh)
}
DATA<-BSdata;gm1<-BSgm

cl.cores <- detectCores()
numWorkers=floor(cl.cores/8)
cl <- makeCluster(numWorkers)
registerDoParallel(cl)
res <- foreach(iter=1:S,.packages=c('foreach','parallel')) %dopar% calculate_gvcm_parallel(iter)
stopCluster(cl)
stopImplicitCluster()

reN=nrow(res[[1]])
boo_result_K<-array(NA,dim=c(S,reN,P+3))
for(i in 1:S){boo_result_K[i,,]<-res[[i]]}
dput(boo_result_K,"/home/csr5/gvcm phone/result/bootstrap/K 55 d9reorder Norm1.5")
# boos<-dget("/home/csr5/gvcm phone/result/bootstrap/K 45 1-200")

source('CauculateBootstrap.R')
dput(Kre,"~/Documents/gvcm phone/result/plot/K estimater d09 norm1 75")
dput(Kre_dp,"~/Documents/gvcm phone/result/plot/Kdp estimater d09 norm1 75")