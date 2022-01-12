library(dplyr)
library(xtable)
library(parallel)
library(ggplot2)
library(flexmix)
source("../EM_functions_MAR.R")
set.seed(1222)


invlogit = function(x) exp(x)/(1+exp(x))

maxT = 8
n = c(100,200,300,400,500)

its = 1000

t3=0.01
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)

t1_hat_ME = list()

t1_hat_CP = list()

t2=invlogit(b0+b2*(0:(maxT-1)))


for(j in 1:length(n)){

  t1_hat_ME[[j]] = rep(NA,its)
  
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
    
    Cov <- Cov[c(3,4,1,2),c(3,4,1,2)]
    
    b1 <- coefs[4]
    
    dq = c(0,invlogit(b1)*(1-invlogit(b1)),0,0)
    ME1 = qnorm(.975)*sqrt(Cov[2,2])
    t1_hat_ME[[j]][it] = qnorm(.975)*sqrt(dq%*%Cov%*%dq)
    t1_hat_CP[[j]][it] = as.numeric(logit(t1)>b1-ME1 & logit(t1)<b1+ME1)
    
  }
      
  df_t1_me= data.frame(t1_hat=t1_hat_ME[[j]])
  write.csv(df_t1_me,paste("SS/df_t1_ME_fm",n[j],".csv",sep=""))
    
  df_t1_cp= data.frame(t1_hat=t1_hat_CP[[j]])
    
  write.csv(df_t1_cp,paste("SS/df_t1_cp_fm",n[j],".csv",sep=""))
    

}
  
  
