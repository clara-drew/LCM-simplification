
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

signif = function(boot_ci,i){
  return(1-coverp(0,boot_ci,i))
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
  ids=sample(1:n[j],n[j],replace=T)
  dfB=data.frame(t=NULL,ids=NULL,y=NULL)
  for(i in 1:n[j]){  
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df[df$ids==ids[i],])
    ni = dim(df[df$ids==ids[i],])[1]
    dfB$ids[(nl+1):(nl+ni)]=i
  }
  
  z_samp = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>0,1,0)
  if(mean(z_samp)==0){

    return(NA)
  }  
  else{ 
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

    m = glm(t2tilp~c(1:maxT),family=binomial,weights=rep(sum(z_samp),8))
    
    return(m$coefficients)
  }
}


bootPow = function(its=1000,cores=3){
  ests = mclapply(1:its,FUN=iter,mc.cores=cores)
  ests = ests[!is.na(ests)]
  mats = lapply(ests,FUN=function(x) (x-m1$coefficients)%*%t(x-m1$coefficients))
  Var = Reduce("+",mats)/length(ests)
  p = length(ests)/its
  se = sqrt(Var[2,2])
  pow = mean(ifelse(pnorm(abs(m1$coefficients[2]/se),lower.tail=FALSE)<0.025,1,0))*p
  df = data.frame(int=unlist(ests)[c(TRUE,FALSE)],slp=unlist(ests)[c(FALSE,TRUE)])
  if(pow==0){write.csv(df,"boot.csv")}
  return(pow)
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
  if(mean(z_samp)==0){
    return(NA)
  }  
  else{                     
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
  
    t2tilp=ifelse(t2tilp<0,0,t2tilp)
    t2tilp=ifelse(t2tilp>1,1,t2tilp)

    m = glm(t2tilp~c(1:maxT),family=binomial,weights=rep(sum(z_samp),8))
      
    slope=m$coefficients[2]
  }  
  return(slope)
}

bootPow_s= function(its=1000,cores=3){
  ests = mclapply(1:its,FUN=iter_s,mc.cores=cores)
  ests = ests[!is.na(ests)]
  mats = lapply(ests,FUN=function(x) (x-m2$coefficients)%*%t(x-m2$coefficients))
  Var = Reduce("+",mats)/length(ests)
  p = length(ests)/its
  se = sqrt(Var[2,2])
  pow = mean(ifelse(pnorm(abs(m2$coefficients[2]/se),lower.tail=FALSE)<0.025,1,0))*p
}



invlogit = function(x) exp(x)/(1+exp(x))

maxT = 8
n = c(10,15,20,25,30)

its = 1000

t3=0.01
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)
b2=c(-0.1,-0.3,-0.5)

for(k in 1:3){

a_til_pow = list()
a_til_p_pow = list()
a_til_s_pow = list()
a_til_s_p_pow = list()
a_hat_pow = list()
power = list()

a_til = list()
a_til_p = list()
a_til_s = list()
a_til_s_p = list()
a_hat = list()

'

t2=invlogit(b0+b2[k]*(0:(maxT-1)))

start = Sys.time()
 
for(j in 1:length(n)){

  a_til_pow[[j]] = rep(NA,its)
  a_til_s_pow[[j]] = rep(NA,its)
  a_til_p_pow[[j]] = rep(NA,its)
  a_til_s_p_pow[[j]] = rep(NA,its)
  a_hat_pow[[j]] = rep(NA,its)
  power[[j]] = rep(NA,its)

  a_til[[j]] = as.data.frame(matrix(ncol=2,nrow=its,dimnames=list((1:its),c("til_int","til_slp"))),row.names=NULL)
  a_til_s[[j]] = as.data.frame(matrix(ncol=2,nrow=its,dimnames=list((1:its),c("til_s_int","til_s_slp"))),row.names=NULL)
  a_til_p[[j]] = as.data.frame(matrix(ncol=2,nrow=its,dimnames=list((1:its),c("til_p_int","til_p_slp"))),row.names=NULL)
  a_til_s_p[[j]] = as.data.frame(matrix(ncol=2,nrow=its,dimnames=list((1:its),c("til_s_p_int","til_s_p_slp"))),row.names=NULL)
  a_hat[[j]] = as.data.frame(matrix(ncol=2,nrow=its,dimnames=list((1:its),c("hat_int","hat_slp"))),row.names=NULL)
  

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

    m = tryCatch(glm(y[rep(zeta,each=maxT)==1]~t[rep(zeta,each=maxT)==1],family="binomial"),error=function(x) NA)
    if(!is.na(m[1])){
      power[[j]][it]=ifelse(summary(m)$coefficients[2,4]<0.05,1,0)
    }
    else{
      power[[j]][it]=0
    }
       

    # z, simple
    
    df = data.frame(t=t,ids=ids,y=y)
    z_samp = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>0,1,0)
    df_samp = left_join(data.frame(ids=1:n[j],z=z_samp),df)
    
    t1_samp=mean(z_samp)
    if(t1_samp==0){
      a_til[[j]][it,]=c(NA,NA)
      a_til_pow[[j]][it]=0
    }
    else{
      t2_samp = aggregate(y~t,data=df_samp[df_samp$z==1,],FUN=mean)[,2]

      m = glm(t2_samp~c(1:maxT),family=binomial,weights=rep(sum(z_samp),8))

      if(!is.na(m[1])){
        a_til_pow[[j]][it]=ifelse(summary(m)$coefficients[2,4]<0.05,1,0)
        a_til[[j]][it,]=m$coefficients
      }
      else{
        a_til_pow[[j]][it]=0
        a_til[[j]][it,]=c(NA,NA)
      }
      
    }   
    
    # z-star, simple

    z_samp1 = df[!duplicated(df$id),]$y
    z_samp2 = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>1,1,0)
    z_samp_s= ifelse(z_samp1+z_samp2>0,1,0)
    df_samp_s= left_join(data.frame(ids=1:n[j],z=z_samp_s),df)
    t1_samp_s = mean(z_samp_s)
    if(t1_samp_s==0){
      
      a_til_s[[j]][it,]=c(NA,NA)
      a_til_s_pow[[j]][it]=0
    }
    else{
      t2_samp_s= aggregate(df_samp_s[df_samp_s$z==1,]$y,by=list(ids=df_samp_s[df_samp_s$z==1,]$t),FUN=mean)[,2]
      t3_samp_s= mean(df_samp_s[df_samp_s$z==0,]$y)

      m = tryCatch(glm(t2_samp_s~c(1:maxT),family=binomial,weights=rep(sum(z_samp_s),8)),error=function(x) NA)
      if(!is.na(m[1])){
        a_til_s_pow[[j]][it]=ifelse(summary(m)$coefficients[2,4]<0.05,1,0)
        a_til_s[[j]][it,]=m$coefficients
      }
      else{
        a_til_s_pow[[j]][it]=0
        a_til_s[[j]][it,]=c(NA,NA)
      }   
    }
    
    # z, bias corrected

    if(mean(z_samp)==0){
      a_til_p[[j]][it,]=c(NA,NA)
      a_til_p_pow[[j]][it]=0
    }  
    else{
      
      gp = rep(NA,n[j])
  
      for(i in 1:n[j]){
        tps = df_samp[df_samp$ids==i,]$t+1
        ntps = length(tps)
        gp[i]=prod(1-t2_samp[tps],na.rm=T)  
      }
  
      t1tilp=(t1_samp)/(mean(1-gp))
      t2tilp=(t2_samp*sum((1-gp)*t1tilp))/(t1tilp*n[j])

      m1 = tryCatch(glm(t2tilp~c(1:maxT),family=binomial,weights=rep(sum(z_samp),8)),error=function(x) NA)
      if(!is.na(m1[1])){
        a_til_p[[j]][it,]=m1$coefficients
        a_til_p_pow[[j]][it]=ifelse(summary(m1)$coefficients[2,4]<0.05,1,0)
      }
      else{
       a_til_p[[j]][it,]=c(NA,NA)
       a_til_p_pow[[j]][it]=0
     }
    }

    # z-star, bias corrected

    if(mean(z_samp_s)==0){
      
      a_til_s_p[[j]][it,]=c(NA,NA)
      a_til_s_p_pow[[j]][it]=0
    }  
    else{                     
      pz10 = rep(NA,n[j])
      pz11 = rep(NA,n[j])
      pz111 = matrix(nrow=n[j],ncol=maxT)
      pz110 = rep(NA,n[j])
    
      for(i in 1:n[j]){
        tps = df_samp_s[df_samp_s$ids==i,]$t+1
        Ti = length(tps)
        pz10[i] = (1-sum(dbinom(0:1,Ti-1,t3_samp_s)*(1-t3_samp_s)))
        pz11[i] = (1-prod(1-t2_samp_s[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) t2_samp_s[x]*prod(1-t2_samp_s[-x])))))       
        pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-t2_samp_s[tps[-x]]))))
        pz110[i] = 1-(1-t3_samp_s)^(Ti-1)
      }
    
      t1tilp=(t1_samp_s-mean(pz10))/mean(pz11-pz10)
      t2tilp=(t2_samp_s*(mean(pz11)*t1tilp+mean(pz10)*(1-t1tilp))-
        mean(pz110)*t3_samp_s*t1tilp)/(colMeans(pz111,na.rm=T)*t1tilp)
  
      t2tilp=ifelse(t2tilp<0,0,t2tilp)
      t2tilp=ifelse(t2tilp>1,1,t2tilp)

      m2 = tryCatch(glm(t2tilp~c(1:maxT),family=binomial,weights=rep(sum(z_samp_s),8)),error=function(x) NA)
      if(!is.na(m2[1])){
        a_til_s_p[[j]][it,]=m2$coefficients
        a_til_s_p_pow[[j]][it]=ifelse(summary(m2)$coefficients[2,4]<0.05,1,0)
      }
      else{
        a_til_s_p[[j]][it,]=c(NA,NA)
        a_til_s_p_pow[[j]][it]=0
      }   
    }


    # EM algorithm MLE

    if(sum(y)==0){
      a_hat_pow[[j]][it]=0
      a_hat[[j]][it,]=c(NA,NA)
    }
    else{

      EMests = tryCatch(get_EM_est(df),error=function(x) NA)

      if(!is.na(EMests[1])){

        if(invlogit(EMests$b3)<0.00001){
          Cov = tryCatch(calc_Cov_no3(EMests$b0,EMests$b1,EMests$b2,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov, error=function(x) matrix(ncol=3,nrow=3))
        }else{
          Cov = tryCatch(calc_Cov(EMests$b0,EMests$b1,EMests$b2,EMests$b3,EMests$z,df$t,df$y,df$ids,n=n[j])$Cov, error=function(err) matrix(ncol=4,nrow=4))
        }

        a_hat[[j]][it,]=c(EMests[1],EMests[3])
      
        a_hat_pow[[j]][it]=ifelse(EMests$b2+qnorm(.975)*sqrt(Cov[3,3])< 0 | EMests$b2-qnorm(.975)*sqrt(Cov[3,3])>0,1,0)
      }
      
      else{
        a_hat[[j]][it,]=c(NA,NA)
        a_hat_pow[[j]][it]=0
      }  
    }
      
    df_a = cbind(a_til[[j]],a_til_p[[j]],a_til_s[[j]],a_til_p[[j]],a_til_s_p[[j]],a_hat[[j]])
    
    df_a_pow = data.frame(a_til=a_til_pow[[j]],a_til_p=a_til_p_pow[[j]],a_til_s=a_til_s_pow[[j]],a_til_s_p=a_til_s_p_pow[[j]],a_hat=a_hat_pow[[j]])

    df_a_pow[is.na(df_a_pow)]=0

    write.csv(df_a_pow,paste("SS/df_a_pow_",n[j],"_",k,".csv",sep=""))

  }

}

print(Sys.time()-start)


for(i in 1:5){
  a_hat_pow[[i]][is.na(a_hat_pow[[i]])]=0
}

'

for(i in 1:length(n)){
  df = read.csv(paste("SS/df_a_pow_",n[i],"_",k,".csv",sep=""))
  a_til_pow[[i]]=df$a_til
  a_til_p_pow[[i]]=df$a_til_p
  a_til_s_pow[[i]]=df$a_til_s
  a_til_s_p_pow[[i]]=df$a_til_s_p
  a_hat_pow[[i]]=df$a_hat 
}
 

# Power simulation

df = data.frame(power=c(unlist(lapply(a_til_pow,FUN=mean)),unlist(lapply(a_til_s_pow,FUN=mean)),unlist(lapply(a_til_p_pow,FUN=mean)),unlist(lapply(a_til_s_p_pow,FUN=mean)),unlist(lapply(a_hat_pow,FUN=mean))),estimate=c(rep(c(rep("z",5),rep("z*",5)),2),rep("MLE",5)),BC=c(rep("No",10),rep("Yes",10),rep("No",5)),SS=rep(n,5))

power_fun = function(n){  
  X = cbind(rep(1,maxT),1:maxT)
  W = diag(exp(b0+b2[k]*(1:maxT))/(1+exp(b0+b2[k]*(1:maxT)))^2)
  se = sqrt(solve(t(X) %*% W %*% X)[2,2])
  return(sum(pnorm(abs(b2[k])/(se/sqrt(1:n))-qnorm(.975))*dbinom(1:n,n,t1)))
}

df2 = data.frame(power=unlist(lapply(n,power_fun)),SS=n)


p=ggplot(df,aes(x=SS,y=power,color=estimate,linetype=BC)) + geom_line() + ylab("Power") + xlab("Sample Size") + labs(color="Estimate",linetype="Bias-Corrected") + geom_line(data=df2,mapping=aes(x=SS,y=power),color="Black",linetype="dashed")

ggsave(paste("slope_power_",k,".png",sep=""),p)

}
