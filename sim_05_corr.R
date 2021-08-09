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
  ids=sample(1:n,n,replace=T)
  dfB=data.frame(t=NULL,ids=NULL,y=NULL)
  for(i in 1:n){
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df[df$ids==ids[i],])
    ni = dim(df[df$ids==ids[i],])[1]
    dfB$ids[(nl+1):(nl+ni)]=i
  }
  
  z_samp = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>0,1,0)
  ids = 1:100
  df_boot = left_join(data.frame(ids=ids,z=z_samp),dfB)
  
  gp = rep(NA,length(ids))
  
  t1_samp = mean(z_samp)
  t2_samp = aggregate(y~t,data=df_boot[df_boot$z==1,],FUN=mean)[,2]
  
  for(i in 1:n){
    tps = df_boot[df_boot$ids==i,]$t+1
    ntps = length(tps)
    gp[i]=prod(1-t2_samp[tps],na.rm=T)  
  }
  
  t1tilp=(t1_samp)/(mean(1-gp))
  t2tilp=(t2_samp*sum((1-gp)*t1tilp))/(t1tilp*n)
  
  return(list(t1tilp=t1tilp,t2tilp=t2tilp))
}

bootCI = function(its=1000,cores=3){
  est_list = mclapply(1:its,FUN=iter,mc.cores=cores)
  N = length(unlist(est_list))/its
  CIs = list()
  ind = rep(FALSE,N)
  for(i in 1:N){
    indi = ind
    indi[i] = TRUE
    CIs[[i]] = quantile(unlist(est_list)[indi],probs=c(0.025,0.975),na.rm=T)
  }
  
  return(CIs)
}

iter_s= function(x){

  ids=sample(1:n,n,replace=T)
  dfB=data.frame(t=NULL,ids=NULL,y=NULL)
  for(i in 1:n){
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df[df$ids==ids[i],])
    ni = dim(df[df$ids==ids[i],])[1]
    dfB$ids[(nl+1):(nl+ni)]=i
  }
  

  z_samp1 = dfB[!duplicated(dfB$id),]$y
  z_samp2 = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>1,1,0)
  z_samp = ifelse(z_samp1+z_samp2>0,1,0)
  df_boot = left_join(data.frame(ids=1:n,z=z_samp),dfB)
  t1_samp=mean(z_samp)
  t2_samp = aggregate(df_boot[df_boot$z==1,]$y,by=list(ids=df_boot[df_boot$z==1,]$t),FUN=mean)[,2]
  t3_samp = mean(df_boot[df_boot$z==0,]$y)

  pz10 = rep(NA,n)
  pz11 = rep(NA,n)
  pz111 = matrix(nrow=n,ncol=maxT[j])
  pz110 = rep(NA,n)
    
  for(i in 1:n){
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
  t3tilp=(t3_samp*(mean(1-pz11)*t1tilp+mean(1-pz10)*(1-t1tilp))-
        mean(colMeans(1-pz111,na.rm=T)*t2tilp*(1-t1tilp)))/mean((1-pz110)*t1tilp)
  
 
  return(list(t1tilp=t1tilp,t2tilp=t2tilp,t3tilp=t3tilp))
}

bootCI_s= function(its=1000,cores=3){
  est_list = mclapply(1:its,FUN=iter_s,mc.cores=cores)
  N = length(unlist(est_list))/its
  CIs = list()
  ind = rep(FALSE,N)
  for(i in 1:N){
    indi = ind
    indi[i] = TRUE
    CIs[[i]] = quantile(unlist(est_list)[indi],probs=c(0.025,0.975),na.rm=T)
  }
  
  return(CIs)
}



invlogit = function(x) exp(x)/(1+exp(x))

maxT = c(3)
n = 100

its = 1000

t3=0.05
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)
sigma=1

t1_til = list()
t1_til_p = list()
t1_til_s = list()
t1_til_s_p = list()
t1_hat = list()

t2_til = list()
t2_til_p = list()
t2_til_s = list()
t2_til_s_p = list()
t2_hat = list()

t3_til_s= list()
t3_til_s_p = list()
t3_hat = list()

t1_til_CP = list()
t1_til_p_CP = list()
t1_til_s_CP = list()
t1_til_s_p_CP = list()
t1_hat_CP = list()

t2_til_CP = list()
t2_til_p_CP = list()
t2_til_s_CP = list()
t2_til_s_p_CP = list()
t2_hat_CP = list()

t3_til_s_CP= list()
t3_til_s_p_CP= list()
t3_hat_CP= list()

start = Sys.time()

  
for(j in 1:length(maxT)){

  t2=invlogit(b0+b2[j]*(0:(maxT[j]-1)))

  t1_til[[j]] = rep(NA,its)
  t1_til_s[[j]] = rep(NA,its)
  t2_til[[j]] = matrix(ncol=maxT[j],nrow=its)
  t2_til_s[[j]] = matrix(ncol=maxT[j],nrow=its)
  t1_til_p[[j]] = rep(NA,its)
  t1_til_s_p[[j]] = rep(NA,its)
  t2_til_p[[j]] = matrix(ncol=maxT[j],nrow=its)
  t2_til_s_p[[j]] = matrix(ncol=maxT[j],nrow=its)
  t1_hat[[j]] = rep(NA,its)
  t2_hat[[j]] = matrix(ncol=maxT[j],nrow=its)
  
  t3_hat[[j]] = rep(NA,its)
  t3_til_s[[j]] = rep(NA,its)
  t3_til_s_p[[j]] = rep(NA,its)
    
  t1_til_CP[[j]] = rep(NA,its)
  t1_til_p_CP[[j]] = rep(NA,its)
  t1_til_s_CP[[j]] = rep(NA,its)
  t1_til_s_p_CP[[j]] = rep(NA,its)
  t2_til_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  t2_til_p_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  t2_til_s_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  t2_til_s_p_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  t1_hat_CP[[j]] = rep(NA,its)
  t2_hat_CP[[j]] = matrix(ncol=maxT[j],nrow=its)

  t3_hat_CP[[j]] = rep(NA,its)
  t3_til_s_CP[[j]] = rep(NA,its)
  t3_til_s_p_CP[[j]] = rep(NA,its)
    
  for(it in 1:its){
    
    rand=rnorm(n,0,sigma)

    t2i=matrix(nrow=n,ncol=maxT[j])
    for(i in 1:n){
      t2i[i,]=invlogit(b0+b2[j]*(0:(maxT[j]-1))+rand[i])
    }

    t3i=invlogit(logit(t3)+rand)
    
    zeta=rbinom(n,1,t1)
    t=rep(0:(maxT[j]-1),n)
    ids=rep(1:n,each=maxT[j])
    y=rep(NA,maxT[j]*n)
    for(i in 1:n){
      for(l in 0:(maxT[j]-1)){
        if(zeta[i]==1){y[ids==i & t==l]=rbinom(1,1,t2i[i,l+1])}
        else{y[ids==i & t==l]=rbinom(1,1,t3i[i])}
      }
    }

    # z, simple
    
    df = data.frame(t=t,ids=ids,y=y)
    z_samp = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>0,1,0)
    df_samp = left_join(data.frame(ids=1:n,z=z_samp),df)
   
    t1_samp=mean(z_samp)
    t2_samp = aggregate(y~t,data=df_samp[df_samp$z==1,],FUN=mean)[,2]    
    t1_til[[j]][it]=t1_samp
    t2_til[[j]][it,]=t2_samp
    
    t1_til_CP[[j]][it]=cover(t1,t1_samp,n)
    n1=sum(z_samp)
    t2_til_CP[[j]][it,]=cover(t2,t2_samp,n1)

    # z-star, simple

    z_samp1 = df[!duplicated(df$id),]$y
    z_samp2 = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>1,1,0)
    z_samp_s= ifelse(z_samp1+z_samp2>0,1,0)
    df_samp_s= left_join(data.frame(ids=1:n,z=z_samp_s),df)
    t1_samp_s = mean(z_samp_s)
    t2_samp_s = aggregate(y~t,data=df_samp_s[df_samp_s$z==1,],FUN=mean)[,2]
    t3_samp_s = mean(df_samp_s[df_samp_s$z==0,]$y)
    t1_til_s[[j]][it]=t1_samp_s
    t2_til_s[[j]][it,]=t2_samp_s
    t3_til_s[[j]][it]=t3_samp_s

    t1_til_s_CP[[j]][it]=cover(t1,t1_samp_s,n)
    n1=sum(z_samp_s)
    n0=dim(df_samp_s[df_samp_s$z==0,])[1]
    t2_til_s_CP[[j]][it,]=cover(t2,t2_samp_s,n1)
    t3_til_s_CP[[j]][it]=cover(t3,t3_samp_s,n0)

    # z, bias corrected

    gp = rep(NA,n)
    
    for(i in 1:n){
      tps = df_samp[df_samp$ids==i,]$t+1
      gp[i]=prod(1-t2_samp[tps],na.rm=T)  
    }
    
    t1tilp=(t1_samp)/(mean(1-gp))
    t1_til_p[[j]][it]=t1tilp
    t2_til_p[[j]][it,]=(t2_samp*mean(1-gp))
      
    boot_ci = bootCI(its=1000,cores=23)
    t1_til_p_CP[[j]][it]=coverp(t1,boot_ci,1)
    t2_til_p_CP[[j]][it,]=coverp(t2,boot_ci,2:(1+maxT[j]))

    # z-star, bias corrected

    pz10 = rep(NA,n)
    pz11 = rep(NA,n)
    pz111 = matrix(nrow=n,ncol=maxT[j])
    pz110 = rep(NA,n)

    for(i in 1:n){
      tps = df_samp_s[df_samp_s$ids==i,]$t+1
      Ti = length(tps)
      pz10[i] = 1-sum(dbinom(0:1,Ti-1,t3_samp_s)*(1-t3_samp_s))
      pz11[i] = 1-prod(1-t2_samp_s[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) t2_samp_s[x]*prod(1-t2_samp_s[-x]))))      
      pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-t2_samp_s[tps[-x]]))))
      pz110[i] = 1-(1-t3_samp_s)^(Ti-1)
    }
    
    t1tilp=(t1_samp_s-mean(pz10))/mean(pz11-pz10)
    t1_til_s_p[[j]][it]=t1tilp
    t2tilp=(t2_samp_s*(mean(pz11)*t1tilp+mean(pz10)*(1-t1tilp))-
        mean(pz110)*t3_samp_s*t1tilp)/(colMeans(pz111,na.rm=T)*t1tilp)
    t2_til_s_p[[j]][it,]=t2tilp
    t3_til_s_p[[j]][it]=(t3_samp_s*(mean(1-pz11)*t1tilp+mean(1-pz10)*(1-t1tilp))-
        mean(colMeans(1-pz111,na.rm=T)*t2tilp*(1-t1tilp)))/mean((1-pz110)*t1tilp)
      
    boot_ci= bootCI_s(its=1000,cores=23)
      
    t1_til_s_p_CP[[j]][it]=coverp(t1,boot_ci,1)
    t2_til_s_p_CP[[j]][it,]=coverp(t2,boot_ci,2:(1+maxT[j]))
    t3_til_s_p_CP[[j]][it]=coverp(t3,boot_ci,2+maxT[j])

    # EM algorithm MLE

    EMests = get_EM_est(df)
    t1_hat[[j]][it] = invlogit(EMests$b1)
    t2_hat[[j]][it,] = invlogit(EMests$b0+EMests$b2*(0:(maxT[j]-1)))
    t3_hat[[j]][it] = invlogit(EMests$b3)


    if(t3_hat[[j]][it]<0.00001){
      Cov = calc_Cov_no3(EMests$b0,EMests$b1,EMests$b2,EMests$z,df$t,df$y,df$ids,n=n)$Cov
      ME1 = qnorm(.975)*sqrt(Cov[2,2])
      t1_hat_CP[[j]][it] = as.numeric(logit(t1)>EMests$b1-ME1 & logit(t1)<EMests$b1+ME1)
      for(t in 1:maxT[j]){
        dqdb = c(1,0,(t-1))
        ME2 = qnorm(.975)*sqrt(dqdb %*% Cov %*% dqdb)
        t2_hat_CP[[j]][it,t] = as.numeric(logit(t2[t])>EMests$b0+EMests$b2*(t-1)-ME2 & logit(t2[t])<EMests$b0+EMests$b2*(t-1)+ME2)
      }
      t3_hat_CP[[j]][it] = NA
    }
    else{
      Cov = tryCatch(calc_Cov(EMests$b0,EMests$b1,EMests$b2,EMests$b3,EMests$z,df$t,df$y,df$ids,n=n)$Cov, error=function(err) matrix(ncol=4,nrow=4))
      ME1 = qnorm(.975)*sqrt(Cov[2,2])
      t1_hat_CP[[j]][it] = as.numeric(logit(t1)>EMests$b1-ME1 & logit(t1)<EMests$b1+ME1)
      for(t in 1:maxT[j]){
        dqdb = c(1,0,(t-1),0)
        ME2 = qnorm(.975)*sqrt(dqdb %*% Cov %*% dqdb)
        t2_hat_CP[[j]][it,t] = as.numeric(logit(t2[t])>EMests$b0+EMests$b2*(t-1)-ME2 & logit(t2[t])<EMests$b0+EMests$b2*(t-1)+ME2)
      }
      ME3 = qnorm(.975)*sqrt(Cov[4,4])
      t3_hat_CP[[j]][it] = as.numeric(logit(t3[j])>EMests$b3-ME3 & logit(t3[j])<EMests$b3+ME3)
    }

    df_t1 = data.frame(t1_til=t1_til[[j]],t1_til_p=t1_til_p[[j]],t1_til_s=t1_til_s[[j]],t1_til_s_p=t1_til_s_p[[j]],t1_hat=t1_hat[[j]])
    df_t2 = cbind(t2_til[[j]],t2_til_p[[j]],t2_til_s[[j]],t2_til_s_p[[j]],t2_hat[[j]])
    colnames(df_t2)=c(paste("t2_til_",1:maxT[j],sep=""),paste("t2_til_p_",1:maxT[j],sep=""),paste("t2_til_s_",1:maxT[j],sep=""),paste("t2_til_s_p_",1:maxT[j],sep=""),paste("t2_hat_",1:maxT[j],sep=""))
    write.csv(df_t1,paste("ct305/df_","t1_",maxT[j],".csv",sep=""))
    write.csv(df_t2,paste("ct305/df_","t2_",maxT[j],".csv",sep=""))

    df_t1_cp= data.frame(t1_til=t1_til_CP[[j]],t1_til_p=t1_til_p_CP[[j]],t1_til_s=t1_til_s_CP[[j]],t1_til_s_p=t1_til_s_p_CP[[j]],t1_hat=t1_hat_CP[[j]])
    df_t2_cp= cbind(t2_til_CP[[j]],t2_til_p_CP[[j]],t2_til_s_CP[[j]],t2_til_s_p_CP[[j]],t2_hat_CP[[j]])
    colnames(df_t2_cp)=c(paste("t2_til_",1:maxT[j],sep=""),paste("t2_til_p_",1:maxT[j],sep=""),paste("t2_til_s_",1:maxT[j],sep=""),paste("t2_til_s_p_",1:maxT[j],sep=""),paste("t2_hat_",1:maxT[j],sep=""))
    write.csv(df_t1_cp,paste("ct305/df_","t1_cp_",maxT[j],".csv",sep=""))
    write.csv(df_t2_cp,paste("ct305/df_","t2_cp_",maxT[j],".csv",sep=""))    
 
  }


}


