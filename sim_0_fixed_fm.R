library(dplyr)
library(xtable)
library(parallel)
library(ggplot2)
library(flexmix)
source("../EM_functions_MAR.R")
set.seed(122)

#system("unzip t30.zip")

invlogit = function(x) exp(x)/(1+exp(x))

maxT = 3:8
n = 100

its = 1000

t3=0
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)

t1_hat = list()

t2_hat = list()

t3_hat = list()

t1_hat_CP = list()

t2_hat_CP = list()

t3_hat_CP= list()

start = Sys.time()


for(j in 1:length(maxT)){
  
  t2=invlogit(b0+b2[j]*(0:(maxT[j]-1)))
  
  t1_hat[[j]] = rep(NA,its)
  t2_hat[[j]] = matrix(ncol=maxT[j],nrow=its)
  
  t3_hat[[j]] = rep(NA,its)
  
  t1_hat_CP[[j]] = rep(NA,its)
  t2_hat_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  
  t3_hat_CP[[j]] = rep(NA,its)

  for(it in 1:its){
    zeta=rbinom(n,1,t1)
    t=rep(0:(maxT[j]-1),n)
    ids=rep(1:n,each=maxT[j])
    y=rep(NA,maxT[j]*n)
    for(i in 1:n){
      for(l in 0:(maxT[j]-1)){
        if(zeta[i]==1){y[ids==i & t==l]=rbinom(1,1,t2[l+1])}
        else{y[ids==i & t==l]=rbinom(1,1,t3)}
      }
    }
    
    # EM algorithm MLE
    
    df <- data.frame(y,t,ids)
    
    df_short <- df %>% group_by(ids) %>% filter(row_number() == 1)
    df$Ez <- NA
    
    k=1
    for(i in 1:dim(df_short)[1]){
      times<-as.numeric(table(df$ids)[i])
      df$Ez[k:(k+times-1)] <- rep(df_short$y[i]+1,times)
      k<-k+times
    }
    
    m1 <- flexmix(cbind(y,1-y) ~ 1 | ids, data=df, k=2, cluster=df$Ez,
                  model = FLXMRglmfix(nested = list(k=c(1,1), formula = c(~0,~t)),family="binomial"))
    
    refit_m1 <- refit(m1)
    
    coefs <- refit_m1@coef
    Cov <- refit_m1@vcov
    SEs <- sqrt(diag(Cov))
    
    t1_hat[[j]][it] = invlogit(coefs[4])
    t2_hat[[j]][it,] = invlogit(coefs[3]+coefs[1]*(0:(maxT[j]-1)))
    t3_hat[[j]][it] = invlogit(coefs[2])
    
    ME1 = qnorm(.975)*SEs[4]
    t1_hat_CP[[j]][it] = as.numeric(logit(t1)>coefs[4]-ME1 & logit(t1)<coefs[4]+ME1)
    for(t in 1:maxT[j]){
      dqdb = c((t-1),0,1,0)
      ME2 = qnorm(.975)*sqrt(dqdb %*% Cov %*% dqdb)
      t2_hat_CP[[j]][it,t] = as.numeric(logit(t2[t])>coefs[3]+coefs[1]*(t-1)-ME2 & logit(t2[t])<coefs[3]+coefs[1]*(t-1)+ME2)
    }
    ME3 = qnorm(.975)*SEs[2]
    t3_hat_CP[[j]][it] = as.numeric(logit(t3[j])>coefs[2]-ME3 & logit(t3[j])<coefs[2]+ME3)
    
    df_t1 = data.frame(t1_hat=t1_hat[[j]])
    df_t2 = cbind(t2_hat[[j]])
    colnames(df_t2)=c(paste("t2_hat_",1:maxT[j],sep=""))
    write.csv(df_t1,paste("t30/df_","t1_",maxT[j],"_fm.csv",sep=""))
    write.csv(df_t2,paste("t30/df_","t2_",maxT[j],"_fm.csv",sep=""))
    df_t1_cp= data.frame(t1_hat=t1_hat_CP[[j]])
    df_t2_cp= cbind(t2_hat_CP[[j]])
    colnames(df_t2_cp)=c(paste("t2_hat_",1:maxT[j],sep=""))
    write.csv(df_t1_cp,paste("t30/df_","t1_cp_",maxT[j],"_fm.csv",sep=""))
    write.csv(df_t2_cp,paste("t30/df_","t2_cp_",maxT[j],"_fm.csv",sep=""))
    
  }
  
}


