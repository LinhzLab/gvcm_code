boo0<-boo_result_K; data1<-BSd
reK<-K.result.scaled
idc<-which(boo0[,1,P+3]<10^-4);lidc<-length(idc);P=ncol(data1$X);N=nrow(data1$X)
boo1<-boo0[idc,,]
booR<-array(NA,dim=c(lidc,N,P+2));sde<-matrix(NA,N,P+2)
Gh<-rowSums(data1$X*reK[,1:P])
Vh<-reK[,P+2]

booR_dp<-array(NA,dim=c(lidc,200,P+2));sde_dp<-matrix(NA,200,P+2);reK_dp<-matrix(NA,200,P+2)
boo_dp1<-seq(range(data1$U)[1],range(data1$U)[2],diff(range(data1$U))/199)
boo_dp2<-seq(range(Gh)[1],range(Gh)[2],diff(range(Gh))/199)
boo_dp3<-seq(range(reK[,P+1])[1],range(reK[,P+1])[2],diff(range(reK[,P+1]))/199)
for(p in 1:P){reK_dp[,p]<-approx(data1$U,reK[,p],xout=boo_dp1)$y}
reK_dp[,P+1]<-approx(Gh,reK[,P+1],xout=boo_dp2)$y
reK_dp[,P+2]<-approx(reK[,P+1],reK[,P+2],xout=boo_dp3)$y


for(i in 1:lidc)
{
 for(j in 1:P)
 {
  booR[i,,j]<-approx(DATA[[idc[i]]]$U,boo1[i,,j],xout=data1$U)$y
  booR_dp[i,,j]<-approx(DATA[[idc[i]]]$U,boo1[i,,j],xout=boo_dp1)$y
 }
  booR[i,,P+1]<-approx(rowSums(DATA[[idc[i]]]$X*boo1[i,,1:P]),boo1[i,,P+1],xout=Gh)$y
  booR_dp[i,,P+1]<-approx(rowSums(DATA[[idc[i]]]$X*boo1[i,,1:P]),boo1[i,,P+1],xout=boo_dp2)$y
  booR[i,,P+2]<-approx(boo1[i,,P+1],boo1[i,,P+2],xout=reK[,P+1])$y
  booR_dp[i,,P+2]<-approx(boo1[i,,P+1],boo1[i,,P+2],xout=boo_dp3)$y
  
}

for(p in 1:(P+2))##
{
  for(i in 1:N){ sde[i,p]<-sd(booR[,i,p],na.rm=T) }
  for(i in 1:200){ sde_dp[i,p]<-sd(booR_dp[,i,p],na.rm=T) }
}

upband<-(reK[,1:(P+2)]+1.96*sde)
loband<-(reK[,1:(P+2)]-1.96*sde)
Kre<-list(reK,Gh,upband,loband)

upband_dp<-(reK_dp[,1:(P+2)]+1.96*sde_dp)
loband_dp<-(reK_dp[,1:(P+2)]-1.96*sde_dp)
Kre_dp<-list(reK_dp,boo_dp2,upband_dp,loband_dp)