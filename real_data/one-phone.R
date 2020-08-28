rm(list = ls())
dev.off()
setwd('/home/csr5/gvcm phone')
source('preparefunction.R')
KR<-dget("F:/ÂÛÎÄ/gvcm/´úÂë/gvcmphone/plot/K estimater d09reorder norm1 55");data1<-d9reorder
P=ncol(data1$X);N=nrow(data1$X)
reK<-KR[[1]];Gh<-KR[[2]];upband0<-KR[[3]];loband0<-KR[[4]];upband<-upband0;loband<-loband0

Xnames<-c(names(data.frame(d1$X))[-c(1,2)],"MonthlyIncome")
Unames<-names(data.frame(d1$X))[1]

U<-d1$X[,1]
#U<-data1$U

ifelse(nrow(upband)==200,boo_dp1<-seq(range(U)[1],range(U)[2],diff(range(U))/199),boo_dp1<-U)

Gh<-as.numeric(Gh)
change1=which(Gh==max(Gh[Gh<26]))
boo_dp1<-as.numeric(boo_dp1)
change=which(boo_dp1==max(boo_dp1[boo_dp1<1525.525]))[1]
par(mfrow=c(2,3))
#for(p in c(2,1,3:P))
#{
 plot(boo_dp1[order(boo_dp1)],reK[order(boo_dp1),1],ylim=c(min(loband[,1]),max(upband[,1])),type="l",xlab=Unames,ylab="impact",main=list("(a).Age",cex=1.5,font=1))
 lines(boo_dp1[order(boo_dp1)],upband[order(boo_dp1),1],lty=10)
 lines(boo_dp1[order(boo_dp1)],loband[order(boo_dp1),1],lty=10) 
 lines(boo_dp1,rep(0,nrow(upband0)),lty=3)
 lines(rep(boo_dp1[change],2),range(c(min(loband[,1]),max(upband[,1]))),lty=4)
 
 plot(boo_dp1[order(boo_dp1)],reK[order(boo_dp1),2],ylim=c(min(loband[,2]),max(upband[,2])),type="l",xlab=Unames,ylab="impact",main=list("(b).TongDunScore",cex=1.5,font=1))
 lines(boo_dp1[order(boo_dp1)],upband[order(boo_dp1),2],lty=10)
 lines(boo_dp1[order(boo_dp1)],loband[order(boo_dp1),2],lty=10) 
 lines(boo_dp1,rep(0,nrow(upband0)),lty=3)
 lines(rep(boo_dp1[change],2),range(c(min(loband[,2]),max(upband[,2]))),lty=4)
 
 plot(boo_dp1[order(boo_dp1)],reK[order(boo_dp1),3],ylim=c(min(loband[,3]),max(upband[,3])),type="l",xlab=Unames,ylab="impact",main=list("(c).FirstPayRatio",cex=1.5,font=1))
 lines(boo_dp1[order(boo_dp1)],upband[order(boo_dp1),3],lty=10)
 lines(boo_dp1[order(boo_dp1)],loband[order(boo_dp1),3],lty=10) 
 lines(boo_dp1,rep(0,nrow(upband0)),lty=3)
 lines(rep(boo_dp1[change],2),range(c(min(loband[,3]),max(upband[,3]))),lty=4)
 
 plot(boo_dp1[order(boo_dp1)],reK[order(boo_dp1),4],ylim=c(min(loband[,4]),max(upband[,4])),type="l",xlab=Unames,ylab="impact",main=list("(d).PlantformCount",cex=1.5,font=1))
 lines(boo_dp1[order(boo_dp1)],upband[order(boo_dp1),4],lty=10)
 lines(boo_dp1[order(boo_dp1)],loband[order(boo_dp1),4],lty=10) 
 lines(boo_dp1,rep(0,nrow(upband0)),lty=3)
 lines(rep(boo_dp1[change],2),range(c(min(loband[,4]),max(upband[,4]))),lty=4)
 
 plot(boo_dp1[order(boo_dp1)],reK[order(boo_dp1),5],ylim=c(min(loband[,5]),max(upband[,5])),type="l",xlab=Unames,ylab="impact",main=list("(e).MonthlyIncome",cex=1.5,font=1))
 lines(boo_dp1[order(boo_dp1)],upband[order(boo_dp1),5],lty=10)
 lines(boo_dp1[order(boo_dp1)],loband[order(boo_dp1),5],lty=10) 
 lines(boo_dp1,rep(0,nrow(upband0)),lty=3)
 lines(rep(boo_dp1[change],2),range(c(min(loband[,5]),max(upband[,5]))),lty=4)
#}
plot(Gh[order(Gh)],reK[order(Gh),P+1],ylim=c(min(loband[,P+1]),max(upband[,P+1])),type="l",xlab="x",ylab="g(x)",main=list("(f).link function",cex=1.5,font=1))
lines(Gh[order(Gh)],upband[order(Gh),P+1],lty=2)
lines(Gh[order(Gh)],loband[order(Gh),P+1],lty=2)
lines(rep(Gh[change1],2),range(c(min(loband[,P+1]),max(upband[,P+1]))),lty=4)
#par(mfrow=c(3,1))
#par(mfg=c(3,1,3,1))
plot(reK[order(reK[,P+1]),P+1],reK[order(reK[,P+1]),P+2],ylim=range(c(loband[,P+2],upband[,P+2],reK[,P+1]*(1-reK[,P+1]))),type="l",xlab="x",ylab="v(x)",main=list("(g).variance function",cex=1.5,font=1))
lines(reK[order(reK[,P+1]),P+1],upband[order(reK[,P+1]),P+2],lty=2)
lines(reK[order(reK[,P+1]),P+1],loband[order(reK[,P+1]),P+2],lty=2)
lines(reK[order(reK[,P+1]),P+1],(reK[,P+1]*(1-reK[,P+1]))[order(reK[,P+1])],lty=2,col="red")
#gvcm
par(mfrow=c(2,3)); 
gm0<-dget("initial estimater d9reorder");P=ncol(gm0[,1:(ncol(gm0)-2)]);N=nrow(gm0)
Unames<-names(data.frame(d1$X))[1]
plot(boo_dp1[order(boo_dp1)],gm0[order(d9reorder$U),1],type="l",xlab=Unames,ylab="impact",main=list("(a).Age",cex=1.5,font=1))
plot(boo_dp1[order(boo_dp1)],gm0[order(d9reorder$U),2],type="l",xlab=Unames,ylab="impact",main=list("(b).TongDunScore",cex=1.5,font=1))
plot(boo_dp1[order(boo_dp1)],gm0[order(d9reorder$U),3],type="l",xlab=Unames,ylab="impact",main=list("(c).FirstPayRatio",cex=1.5,font=1))
plot(boo_dp1[order(boo_dp1)],gm0[order(d9reorder$U),4],type="l",xlab=Unames,ylab="impact",main=list("(d).PlantformCount",cex=1.5,font=1))
plot(boo_dp1[order(boo_dp1)],gm0[order(d9reorder$U),5],type="l",xlab=Unames,ylab="impact",main=list("(e).MonthlyIncome",cex=1.5,font=1))
plot(rowSums(gm0[,1:P]*d9reorder$X)[order(rowSums(gm0[,1:P]*d9reorder$X))],gm0[order(rowSums(gm0[,1:P]*d9reorder$X)),1+P],type="l",xlab="x",ylab="g(x)",main=list("(f).link function",cex=1.5,font=1))

savePlot("gvcm.pdf",type = "pdf")


