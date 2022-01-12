library(dplyr)
library(xtable)
library(parallel)
library(ggplot2)
source("../EM_functions_MAR.R")
set.seed(1222)

cover = function(p,phat,n){
  me = sqrt(phat*(1-phat)/n)
  return(ifelse(p>phat-qnorm(.975)*me & p<phat+qnorm(.975)*me,1,0))
}

coverp = function(p,boot_ci,i){
  cover=c()
  for(k in 1:length(i)){
    cover=c(cover,ifelse(p[k]>boot_ci[[i[k]]][1] & p[k]<boot_ci[[i[k]]][2],1,0))
  }
  return(cover)
}

signif = function(boot_ci,i,j){
  return(ifelse(coverp(boot_ci[[j]][1],boot_ci,i)+coverp(boot_ci[[j]][2],boot_ci,i)<1,1,0))
}

signif_EM = function(EMests,Cov,t){
  t2=invlogit(EMests$b0+EMests$b2*(1:maxT))
  s=c(1,t)
  se=rep(NA,2)
  for(i in 1:2){
    dqdb0 = invlogit(EMests$b0+s[i]*EMests$b2)*(1-invlogit(EMests$b0+s[i]*EMests$b2))
    if(dim(Cov)[1]==4){
      dqdbt = c(dqdb0,0,s[i]*dqdb0,0)
    }
    else{
      dqdbt = c(dqdb0,0,s[i]*dqdb0)
    }
    se[i]=sqrt(dqdbt %*% Cov %*% dqdbt)
  }
  lower1=t2[1]-qnorm(.975)*se[1]
  lowert=t2[t]-qnorm(.975)*se[2]
  upper1=t2[1]+qnorm(.975)*se[1]
  uppert=t2[t]+qnorm(.975)*se[2]
  return(ifelse((lower1>lowert & lower1<uppert) | (upper1>lowert & upper1<uppert),0,1))
         
}

iter = function(x){
  all0=TRUE
  while(all0==TRUE){
    ids=sample(1:n[j],n[j],replace=T)
    dfB=data.frame(t=NULL,ids=NULL,y=NULL)
    for(i in 1:n[j]){  
      nl = dim(dfB)[1]
      dfB = rbind(dfB,df[df$ids==ids[i],])
      ni = dim(df[df$ids==ids[i],])[1]
      dfB$ids[(nl+1):(nl+ni)]=i
    }
  
    z_samp = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>0,1,0)
    if(mean(z_samp)!=0) all0=FALSE
  }
  ids = 1:n[j]
  df_boot = left_join(data.frame(ids=ids,z=z_samp),dfB)
  
  gp = rep(NA,length(ids))
  
  t1_samp = mean(z_samp)
  t2_samp = aggregate(y~t,data=df_boot[df_boot$z==1,],FUN=mean)[,2]
  
  for(i in 1:n[j]){
    tps = df_boot[df_boot$ids==i,]$t+1
    ntps = length(tps)
    gp[i]=prod(1-t2_samp[tps],na.rm=T)  
  }
  
  t1tilp=(t1_samp)/(mean(1-gp))
  t2tilp=(t2_samp*sum((1-gp)*t1tilp))/(t1tilp*n[j])
  
  return(list(t1tilp=t1tilp,t2tilp=t2tilp))
}

bootCI = function(its=1000,cores=3){
  est_list = mclapply(1:its,FUN=iter,mc.cores=cores)
  N = length(unlist(est_list))/its
  CIs = list()
  MEs = list()
  ind = rep(FALSE,N)
  for(i in 1:N){
    indi = ind
    indi[i] = TRUE
    CIs[[i]] = quantile(unlist(est_list)[indi],probs=c(0.025,0.975),na.rm=T)
    MEs[[i]] = qnorm(.975)*sd(unlist(est_list)[indi])
  }
 
  return(list(CIs=CIs,MEs=MEs))
}


iter_s= function(x){
  all0=TRUE
  while(all0==TRUE){
    ids=sample(1:n[j],n[j],replace=T)
    dfB=data.frame(t=NULL,ids=NULL,y=NULL)
    for(i in 1:n[j]){
      nl = dim(dfB)[1]
      dfB = rbind(dfB,df[df$ids==ids[i],])
      ni = dim(df[df$ids==ids[i],])[1]
      dfB$ids[(nl+1):(nl+ni)]=i
    }
  

    z_samp1 = dfB[!duplicated(dfB$id),]$y
    z_samp2 = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>1,1,0)
    z_samp = ifelse(z_samp1+z_samp2>0,1,0)
    if(mean(z_samp)!=0) all0=FALSE
  }  
  df_boot = left_join(data.frame(ids=1:n[j],z=z_samp),dfB)
  t1_samp=mean(z_samp)
  t2_samp = aggregate(df_boot[df_boot$z==1,]$y,by=list(ids=df_boot[df_boot$z==1,]$t),FUN=mean)[,2]
  t3_samp = mean(df_boot[df_boot$z==0,]$y)

  pz10 = rep(NA,n[j])
  pz11 = rep(NA,n[j])
  pz111 = matrix(nrow=n[j],ncol=maxT)
  pz110 = rep(NA,n[j])
    
  for(i in 1:n[j]){
    tps = df_boot[df_boot$ids==i,]$t+1
    Ti = length(tps)
    pz10[i] = (1-sum(dbinom(0:1,Ti-1,t3_samp)*(1-t3_samp)))
    pz11[i] = (1-prod(1-t2_samp[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) t2_samp[x]*prod(1-t2_samp[-x])))))       
    pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-t2_samp[tps[-x]]))))
        pz110[i] = 1-(1-t3_samp)^(Ti-1)
  }
    
  t1tilp=(t1_samp-mean(pz10))/mean(pz11-pz10)
  t2tilp=(t2_samp*(mean(pz11)*t1tilp+mean(pz10)*(1-t1tilp))-
        mean(pz110)*t3_samp*t1tilp)/(colMeans(pz111,na.rm=T)*t1tilp)
  
 
  return(list(t1tilp=t1tilp,t2tilp=t2tilp))
}

bootCI_s= function(its=1000,cores=3){
  est_list = mclapply(1:its,FUN=iter_s,mc.cores=cores)
  N = length(unlist(est_list))/its
  CIs = list()
  MEs = list()
  ind = rep(FALSE,N)
  for(i in 1:N){
    indi = ind
    indi[i] = TRUE
    CIs[[i]] = quantile(unlist(est_list)[indi],probs=c(0.025,0.975),na.rm=T)
    MEs[[i]] = sd(unlist(est_list)[indi])
  }
  
  return(list(CIs=CIs,MEs=MEs))
}



invlogit = function(x) exp(x)/(1+exp(x))

maxT = 8
n = c(50,125,200,275,350,425)

its = 1000

t3=0.01
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)

a_til_pow4 = list()
a_til_p_pow4 = list()
a_til_s_pow4 = list()
a_til_s_p_pow4 = list()
a_hat_pow4 = list()

a_til_pow6 = list()
a_til_p_pow6 = list()
a_til_s_pow6 = list()
a_til_s_p_pow6 = list()
a_hat_pow6 = list()

a_til_pow8 = list()
a_til_p_pow8 = list()
a_til_s_pow8 = list()
a_til_s_p_pow8 = list()
a_hat_pow8 = list()


t2=invlogit(b0+b2*(0:(maxT-1)))

start = Sys.time()
 
for(j in 1:length(n)){

  a_til_pow4[[j]] = rep(NA,its)
  a_til_s_pow4[[j]] = rep(NA,its)
  a_til_p_pow4[[j]] = rep(NA,its)
  a_til_s_p_pow4[[j]] = rep(NA,its)
  a_hat_pow4[[j]] = rep(NA,its)

  a_til_pow6[[j]] = rep(NA,its)
  a_til_s_pow6[[j]] = rep(NA,its)
  a_til_p_pow6[[j]] = rep(NA,its)
  a_til_s_p_pow6[[j]] = rep(NA,its)
  a_hat_pow6[[j]] = rep(NA,its)

  a_til_pow8[[j]] = rep(NA,its)
  a_til_s_pow8[[j]] = rep(NA,its)
  a_til_p_pow8[[j]] = rep(NA,its)
  a_til_s_p_pow8[[j]] = rep(NA,its)
  a_hat_pow8[[j]] = rep(NA,its)
    

  for(it in 1:its){

    zeta=rbinom(n[j],1,t1)
    t=rep(0:(maxT-1),n[j])
    ids=rep(1:n[j],each=maxT)
    y=rep(NA,maxT*n[j])
    for(i in 1:n[j]){
      for(l in 0:(maxT-1)){
        if(zeta[i]==1){y[ids==i & t==l]=rbinom(1,1,t2[l+1])}
        else{y[ids==i & t==l]=rbinom(1,1,t3)}
      }
    }

    # z, simple
    
    df = data.frame(t=t,ids=ids,y=y)
    write.csv(df,"df_sim.csv")
    z_samp = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>0,1,0)
    df_samp = left_join(data.frame(ids=1:n[j],z=z_samp),df)
    
    t1_samp=mean(z_samp)
    t2_samp = aggregate(y~t,data=df_samp[df_samp$z==1,],FUN=mean)[,2]

    x=as.factor(df_samp[df_samp$z==1 & df_samp$t==0,]$y)
    levels(x)=c(0,1)
    y=as.factor(df_samp[df_samp$z==1 & df_samp$t==3,]$y)
    levels(y)=c(0,1)

    a_til_pow4[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0)

    y=as.factor(df_samp[df_samp$z==1 & df_samp$t==5,]$y)
    levels(y)=c(0,1)

    a_til_pow6[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0)

    y=as.factor(df_samp[df_samp$z==1 & df_samp$t==7,]$y)
    levels(y)=c(0,1)

    a_til_pow8[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0)

    write.csv(df_samp,"sim_data.csv")
    
    # z-star, simple

    z_samp1 = df[!duplicated(df$id),]$y
    z_samp2 = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>1,1,0)
    z_samp_s= ifelse(z_samp1+z_samp2>0,1,0)
    df_samp_s= left_join(data.frame(ids=1:n[j],z=z_samp_s),df)
    t1_samp_s = mean(z_samp_s)
    t2_samp_s= aggregate(df_samp_s[df_samp_s$z==1,]$y,by=list(ids=df_samp_s[df_samp_s$z==1,]$t),FUN=mean)[,2]
    t3_samp_s= mean(df_samp_s[df_samp_s$z==0,]$y)

    x=as.factor(df_samp_s[df_samp_s$z==1 & df_samp_s$t==0,]$y)
    levels(x)=c(0,1)
    y=as.factor(df_samp_s[df_samp_s$z==1 & df_samp_s$t==3,]$y) 
    levels(y)=c(0,1)

    a_til_s_pow4[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0)

    y=as.factor(df_samp_s[df_samp_s$z==1 & df_samp_s$t==5,]$y)
    levels(y)=c(0,1)

    a_til_s_pow6[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0)

    y=as.factor(df_samp_s[df_samp_s$z==1 & df_samp_s$t==7,]$y)
    levels(y)=c(0,1)        

    a_til_s_pow8[[j]][it]=ifelse(mcnemar.test(x,y)$p.value<0.05,1,0) 

    # z, bias corrected

    gp = rep(NA,n[j])
    
    for(i in 1:n[j]){
      tps = df_samp[df_samp$ids==i,]$t+1
      gp[i]=prod(1-t2_samp[tps],na.rm=T)  
    }
    
    t1tilp=(t1_samp)/(mean(1-gp))
    t2tilp=(t2_samp*sum((1-gp)*t1tilp))/(t1tilp*n[j])
      
    boot = bootCI(its=1000,cores=23)
    boot_ci = boot$CIs

    a_til_p_pow4[[j]][it]=signif(boot_ci,2,5)
    a_til_p_pow6[[j]][it]=signif(boot_ci,2,7)
    a_til_p_pow8[[j]][it]=signif(boot_ci,2,9)
 

    # z-star, bias corrected

    pz10 = rep(NA,n[j])
    pz11 = rep(NA,n[j])
    pz111 = matrix(nrow=n[j],ncol=maxT)
    pz110 = rep(NA,n[j])

    for(i in 1:n[j]){
      tps = df_samp_s[df_samp_s$ids==i,]$t+1
      Ti = length(tps)
      pz10[i] = 1-sum(dbinom(0:1,Ti-1,t3_samp_s)*(1-t3_samp_s))
      pz11[i] = 1-prod(1-t2_samp_s[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) t2_samp_s[x]*prod(1-t2_samp_s[-x]))))      
      pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-t2_samp_s[tps[-x]]))))
      pz110[i] = 1-(1-t3_samp_s)^(Ti-1)
    }
    
    t1tilp=(t1_samp_s-mean(pz10))/mean(pz11-pz10)
     
    boot= bootCI_s(its=1000,cores=23)
    boot_ci = boot$CIs

    a_til_s_p_pow4[[j]][it]=signif(boot_ci,2,5)
    a_til_s_p_pow6[[j]][it]=signif(boot_ci,2,7)
    a_til_s_p_pow8[[j]][it]=signif(boot_ci,2,9)

    # EM algorithm MLE

    EMests = get_EM_est(df)

    if(invlogit(EMests$b3)<0.00001){
      Cov = calc_Cov_no3(EMests$b0,EMests$b1,EMests$b2,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov
      a_hat_pow4[[j]][it]=signif_EM(EMests,Cov,4)
      a_hat_pow6[[j]][it]=signif_EM(EMests,Cov,6)
      a_hat_pow8[[j]][it]=signif_EM(EMests,Cov,8)
    }else{
      Cov = tryCatch(calc_Cov(EMests$b0,EMests$b1,EMests$b2,EMests$b3,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov, error=function(err) matrix(ncol=4,nrow=4))
      a_hat_pow4[[j]][it]=signif_EM(EMests,Cov,4)
      a_hat_pow6[[j]][it]=signif_EM(EMests,Cov,6)
      a_hat_pow8[[j]][it]=signif_EM(EMests,Cov,8)
    }
   
    df_a_pow4 = data.frame(a_til=a_til_pow4[[j]],a_til_p=a_til_pow4[[j]],a_til_s=a_til_s_pow4[[j]],a_til_s_p=a_til_s_p_pow4[[j]],a_hat=a_hat_pow4[[j]])

    write.csv(df_a_pow4,paste("SS/df_a_pow4_",n[j],".csv",sep=""))

    df_a_pow6 = data.frame(a_til=a_til_pow6[[j]],a_til_p=a_til_pow6[[j]],a_til_s=a_til_s_pow6[[j]],a_til_s_p=a_til_s_p_pow6[[j]],a_hat=a_hat_pow6[[j]])

    write.csv(df_a_pow6,paste("SS/df_a_pow6_",n[j],".csv",sep=""))

    df_a_pow8 = data.frame(a_til=a_til_pow8[[j]],a_til_p=a_til_pow8[[j]],a_til_s=a_til_s_pow8[[j]],a_til_s_p=a_til_s_p_pow8[[j]],a_hat=a_hat_pow8[[j]])

    write.csv(df_a_pow8,paste("SS/df_a_pow8_",n[j],".csv",sep=""))

  }

}

print(Sys.time()-start)

