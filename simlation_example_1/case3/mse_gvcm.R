library(MASS)
#path <- paste0("/home/LiuJX/gvcmsim/gvcm_example1_case_3.Rdata")
#load(path) 
#########true parameter##########
data_plot <- function(N){
  #########################
  name <-  paste0("N.", N)
  est_g <- List.Kern[[name]]$result_g
  est_bbeta_var <- List.Kern[[name]]$result_bbeta_var
  beta.1.mise <-mean(apply(est_bbeta_var[,,'beta1']- est_bbeta_var[,,'beta1_true'],1,function(x) sum(x^2/length(which(!is.na(x))),na.rm=T)),is.na = T) 
  beta.2.mise <-mean(apply(est_bbeta_var[,,'beta2']- est_bbeta_var[,,'beta2_true'],1,function(x) sum(x^2/length(which(!is.na(x))),na.rm =T)),is.na = T) 
  beta.3.mise <-mean(apply(est_bbeta_var[,,'beta3']- est_bbeta_var[,,'beta3_true'],1,function(x) sum(x^2/length(which(!is.na(x))),na.rm =T)),is.na = T) 
  g.mise <-mean(apply(est_g[,,'g']- est_g[,,'g_true'],1,function(x) sum(x^2/length(which(!is.na(x))),na.rm=T)),is.na = T) 
  Vari.mise <-mean(apply(abs(est_bbeta_var[,,'V']- est_bbeta_var[,,'Variance_true']),1,function(x) sum(x[which(!is.infinite(x))]^2/length(which(!is.infinite(x))),na.rm=T)),is.na = T) 
  total.mise <- sum(beta.1.mise+beta.2.mise+beta.3.mise)
  MISE <-  data.frame(beta.1.mise,beta.2.mise,beta.3.mise,g.mise,Vari.mise,total.mise)
  MISE <- round(MISE,digits=4)
  name_1 <- paste0('gvcm_example1_case',case,'_size',N,'.csv')
  write.csv(MISE,file=name_1,row.names=FALSE)
  ###################plot####################
  set.seed(1)
  U<-runif(N,0,1)
  dp1<-seq(0,1,1/Points)
  #####estimate true beta function#######
  ###fix dense point u############
  fdp<-seq(0,1,1/(N-1));beta.1=fdp^2+1; beta.2=cos(pi*fdp)*cos(pi*fdp)+rep(0.5,N);beta.3=2*sin(pi*fdp)*sin(pi*fdp)-rep(0.5,N)
  bbeta=matrix(cbind(beta.1,beta.2,beta.3),N,P)
  for(i in 1:N){bbeta[i,]=sign(bbeta[i,1])*bbeta[i,]/sqrt(sum(bbeta[i,1:P]^2))}
  beta.true.1=approx(fdp,bbeta[,1],xout=dp1,rule=2:1)$y
  beta.true.2=approx(fdp,bbeta[,2],xout=dp1,rule=2:1)$y
  beta.true.3=approx(fdp,bbeta[,3],xout=dp1,rule=2:1)$y
  x<-mvrnorm(N, rep(0, P), Sigma);bbeta=matrix(cbind(U^2+1,cos(pi*U)*cos(pi*U)+rep(0.5,N),2*sin(pi*U)*sin(pi*U)+rep(0.5,N)),N,P);Xbeta=rowSums(x*bbeta)
  #dp2<-seq(range(Xbeta)[1],range(Xbeta)[2],diff(range(Xbeta))/Points); g.true =linkFun(dp2,case)[[1]]
  dp2<-seq(quantile(Xbeta,probs = 0.01),quantile(Xbeta,probs = 0.99),diff(c(quantile(Xbeta,probs = 0.01),quantile(Xbeta,probs = 0.99)))/Points); g.true<-linkFun(dp2,case)[[1]]
  dp3<-seq(-1,1,2/Points)
  V.true <- rep(0.01,length(dp3))
  #########################
  est_data <- List.Kern[[name]]$result.fix
  beta.1.mean = apply(est_data[,,'beta.1.fix'],2, function(x) mean(x,na.rm= T))
  beta.2.mean = apply(est_data[,,'beta.2.fix'],2, function(x) mean(x,na.rm= T))
  beta.3.mean = apply(est_data[,,'beta.3.fix'],2, function(x) mean(x,na.rm= T))
  beta.1.sd = apply(est_data[,,'beta.1.fix'],2, function(x) sd(x,na.rm= T))
  beta.2.sd = apply(est_data[,,'beta.2.fix'],2, function(x) sd(x,na.rm= T))
  beta.3.sd = apply(est_data[,,'beta.3.fix'],2, function(x) sd(x,na.rm= T))
  g.mean = apply(est_data[,,'g.fix'],2, function(x) mean(x,na.rm= T))
  g.sd = apply(est_data[,,'g.fix'],2, function(x) sd(x,na.rm= T))
  V.mean = apply(est_data[,,'Variance.fix'],2, function(x) mean(x,na.rm= T))
  V.sd = apply(est_data[,,'Variance.fix'],2, function(x) sd(x,na.rm= T))
  
  beta.1.upper = beta.1.mean+1.96*beta.1.sd
  beta.2.upper = beta.2.mean+1.96*beta.2.sd
  beta.3.upper = beta.3.mean+1.96*beta.3.sd
  g.upper = g.mean+1.96*g.sd;
  V.upper = V.mean+1.96*V.sd;
  beta.1.lower = beta.1.mean-1.96*beta.1.sd; beta.2.lower = beta.2.mean-1.96*beta.2.sd; beta.3.lower = beta.3.mean-1.96*beta.3.sd
  g.lower = g.mean-1.96*g.sd;V.lower = V.mean-1.96*V.sd;
  beta1.plot.data <- data.frame(dp1 = dp1, beta.1.mean = beta.1.mean, beta.1.upper = beta.1.upper,beta.1.lower = beta.1.lower,beta.true.1 = beta.true.1)
  beta2.plot.data <- data.frame(dp1 = dp1, beta.2.mean = beta.2.mean, beta.2.upper = beta.2.upper,beta.2.lower = beta.2.lower,beta.true.2= beta.true.2)
  beta3.plot.data <- data.frame(dp1 = dp1, beta.3.mean = beta.3.mean, beta.3.upper = beta.3.upper,beta.3.lower = beta.3.lower,beta.true.3= beta.true.3)
  g.plot.data <- data.frame(dp2 = dp2, g.mean = g.mean, g.upper = g.upper, g.lower = g.lower,g.true = g.true)
  V.plot.data <- data.frame(dp3 = dp3, V.mean = V.mean, V.upper = V.upper, V.lower = V.lower,V.true = V.true)
  
  source('plot_gmcv_1.R',local=TRUE)
  name2 <-  paste0("gvcm_ex1_case_",case,":number_",N,".png")
  ggsave(name2,plot = a,width = 4, height = 16)
  print(MISE)
}
N = 100
while(N < 4001){
  if(N==4000){resulting <- data_plot(N);N = N+5}
  if(N==400){resulting <- data_plot(N);N=N+3600}
  if(N==200){resulting <-data_plot(N);N=N+200}
  if(N==100){resulting <-data_plot(N);N=N+100} 
}
