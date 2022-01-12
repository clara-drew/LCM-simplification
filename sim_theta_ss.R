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

iter = function(x){
  ids=sample(1:n[j],n[j],replace=T)
  dfB=data.frame(t=NULL,ids=NULL,y=NULL)
  for(i in 1:n[j]){
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df[df$ids==ids[i],])
    ni = dim(df[df$ids==ids[i],])[1]
    dfB$ids[(nl+1):(nl+ni)]=i
  }
  
  z_samp = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>0,1,0)
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
  
  return(list(t1tilp=t1tilp))
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
  
 
  return(list(t1tilp=t1tilp))
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
n = c(100,200,300,400,500)

its = 1000

t3=0.01
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)

t1_til_ME = list()
t1_til_p_ME = list()
t1_til_s_ME = list()
t1_til_s_p_ME = list()
t1_hat_ME = list()

t1_til_CP = list()
t1_til_p_CP = list()
t1_til_s_CP = list()
t1_til_s_p_CP = list()
t1_hat_CP = list()

t2=invlogit(b0+b2*(0:(maxT-1)))

start = Sys.time()

  
for(j in 1:length(n)){

  t1_til_ME[[j]] = rep(NA,its)
  t1_til_s_ME[[j]] = rep(NA,its)
  t1_til_p_ME[[j]] = rep(NA,its)
  t1_til_s_p_ME[[j]] = rep(NA,its)
  t1_hat_ME[[j]] = rep(NA,its)
    
  t1_til_CP[[j]] = rep(NA,its)
  t1_til_p_CP[[j]] = rep(NA,its)
  t1_til_s_CP[[j]] = rep(NA,its)
  t1_til_s_p_CP[[j]] = rep(NA,its)
  t1_hat_CP[[j]] = rep(NA,its)
    
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
    z_samp = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>0,1,0)
    df_samp = left_join(data.frame(ids=1:n[j],z=z_samp),df)
    
    t1_samp=mean(z_samp)
    t2_samp = aggregate(y~t,data=df_samp[df_samp$z==1,],FUN=mean)[,2]
    
    t1_til_CP[[j]][it]=cover(t1,t1_samp,n[j])
    t1_til_ME[[j]][it]=qnorm(.975)*sqrt(t1_samp*(1-t1_samp)/n[j])

    # z-star, simple

    z_samp1 = df[!duplicated(df$id),]$y
    z_samp2 = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>1,1,0)
    z_samp_s= ifelse(z_samp1+z_samp2>0,1,0)
    df_samp_s= left_join(data.frame(ids=1:n[j],z=z_samp_s),df)
    t1_samp_s = mean(z_samp_s)
    t2_samp_s= aggregate(df_samp_s[df_samp_s$z==1,]$y,by=list(ids=df_samp_s[df_samp_s$z==1,]$t),FUN=mean)[,2]
    t3_samp_s= mean(df_samp_s[df_samp_s$z==0,]$y)

    t1_til_s_CP[[j]][it]=cover(t1,t1_samp_s,n[j])
    t1_til_s_ME[[j]][it]=qnorm(.975)*sqrt(t1_samp_s*(1-t1_samp_s)/n[j])

    # z, bias corrected

    gp = rep(NA,n[j])
    
    for(i in 1:n[j]){
      tps = df_samp[df_samp$ids==i,]$t+1
      gp[i]=prod(1-t2_samp[tps],na.rm=T)  
    }
    
    t1tilp=(t1_samp)/(mean(1-gp))
      
    boot = bootCI(its=1000,cores=23)
    boot_ci = boot$CIs
    t1_til_p_ME[[j]][it] = unlist(boot$MEs)
    t1_til_p_CP[[j]][it]=coverp(t1,boot_ci,1)

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
    t1_til_s_p_ME[[j]][it]=unlist(boot$MEs)
      
    t1_til_s_p_CP[[j]][it]=coverp(t1,boot_ci,1)

    # EM algorithm MLE

    EMests = get_EM_est(df)

    if(invlogit(EMests$b3)<0.00001){
      Cov = calc_Cov_no3(EMests$b0,EMests$b1,EMests$b2,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov
      dq = c(0,invlogit(EMests$b1)*(1-invlogit(EMests$b1)),0)
      ME1 = qnorm(.975)*sqrt(Cov[2,2])
      t1_hat_ME[[j]][it] = qnorm(.975)*sqrt(dq%*%Cov%*%dq)
      t1_hat_CP[[j]][it] = as.numeric(logit(t1)>EMests$b1-ME1 & logit(t1)<EMests$b1+ME1)
    }
    else{
      Cov = tryCatch(calc_Cov(EMests$b0,EMests$b1,EMests$b2,EMests$b3,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov, error=function(err) matrix(ncol=4,nrow=4))
      dq = c(0,invlogit(EMests$b1)*(1-invlogit(EMests$b1)),0,0)
      ME1 = qnorm(.975)*sqrt(Cov[2,2])
      t1_hat_ME[[j]][it] = qnorm(.975)*sqrt(dq%*%Cov%*%dq)
      t1_hat_CP[[j]][it] = as.numeric(logit(t1)>EMests$b1-ME1 & logit(t1)<EMests$b1+ME1)
    }
   
    df_t1_me= data.frame(t1_til=t1_til_ME[[j]],t1_til_p=t1_til_p_ME[[j]],t1_til_s=t1_til_s_ME[[j]],t1_til_s_p=t1_til_s_p_ME[[j]],t1_hat=t1_hat_ME[[j]])
    write.csv(df_t1_me,paste("SS/df_t1_ME_",n[j],".csv",sep=""))

    df_t1_cp= data.frame(t1_til=t1_til_CP[[j]],t1_til_p=t1_til_p_CP[[j]],t1_til_s=t1_til_s_CP[[j]],t1_til_s_p=t1_til_s_p_CP[[j]],t1_hat=t1_hat_CP[[j]])

    write.csv(df_t1_cp,paste("SS/df_t1_cp_",n[j],".csv",sep=""))

  }
 

}

print(Sys.time()-start)




                                                                                           
  

