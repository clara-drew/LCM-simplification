library(dplyr)
library(xtable)
library(parallel)
library(ggplot2)
library(flexmix)
library(dplyr)
source("../EM_functions_MAR.R")
set.seed(1222)


invlogit = function(x) exp(x)/(1+exp(x))

maxT = 8
n = c(10,15,20,25,30)

its = 1000

t3=0.01
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)
b2=c(-0.1,-0.3,-0.5)

for(l in 1:3){
  
  a_til_pow = list()
  a_til_p_pow = list()
  a_til_s_pow = list()
  a_til_s_p_pow = list()
  a_hat_pow = list()

  t2=invlogit(b0+b2[l]*(0:(maxT-1)))
 
  for(j in 1:length(n)){

    a_hat_pow[[j]] = rep(NA,its)

    for(it in 1:its){
      zeta=rbinom(n[j],1,t1)
      t=rep(0:(maxT-1),n[j])
      ids=rep(1:n[j],each=maxT)
      y=rep(NA,maxT*n[j])
      for(i in 1:n[j]){
        for(m in 0:(maxT-1)){
          if(zeta[i]==1){y[ids==i & t==m]=rbinom(1,1,t2[m+1])}
          else{y[ids==i & t==m]=rbinom(1,1,t3)}
        }
      }

      # EM algorithm MLE

      if(sum(y)==0){
        a_hat_pow[[j]][it]=0
      }
      else{

        df <- data.frame(y,t,ids)
      
        df_short <- df %>% group_by(ids) %>% filter(row_number() == 1)
        df$Ez <- NA
      
        k=1
        for(i in 1:dim(df_short)[1]){
          times<-as.numeric(table(df$ids)[i])
          df$Ez[k:(k+times-1)] <- rep(df_short$y[i]+1,times)
          k<-k+times
        }
      
        m1 <- tryCatch(flexmix(cbind(y,1-y) ~ 1 | ids, data=df, k=2, cluster=df$Ez,
                    model = FLXMRglmfix(nested = list(k=c(1,1), formula = c(~0,~t)),family="binomial")), error=function(x){return(NA)})
      
        if(is.na(m1)){
          a_hat_pow[[j]][it]=0
        }
      
        else{
      
          refit_m1 <- tryCatch(refit(m1), error=function(x){return(NA)})
          
          if(is.na(refit_m1)){
            a_hat_pow[[j]][it]=0
          }
          
          else{
      
            coefs <- refit_m1@coef
            Cov <- refit_m1@vcov
            SEs <- sqrt(diag(Cov))
        
            if(sum(is.nan(SEs))>0){
              a_hat_pow[[j]][it]=0
            }
      
            else{
              a_hat_pow[[j]][it]=ifelse(coefs[1]+qnorm(.975)*SEs[1] < 0 | coefs[1]-qnorm(.975)*SEs[1]>0,1,0)
            }
          }  
        }  
      }
    }  
      
    df_a_pow <- data.frame(a_hat_pow=a_hat_pow[[j]])  

    write.csv(df_a_pow,paste("SS/df_ahat_pow_fm_",n[j],"_",l,".csv",sep=""))

  }



  for(i in 1:5){
    a_hat_pow[[i]][is.na(a_hat_pow[[i]])]=0
  }



  for(i in 1:length(n)){
    df = read.csv(paste("SS/df_a_pow_",n[i],"_",l,".csv",sep=""))
    a_til_pow[[i]]=df$a_til
    a_til_p_pow[[i]]=df$a_til_p
    a_til_s_pow[[i]]=df$a_til_s
    a_til_s_p_pow[[i]]=df$a_til_s_p
  }


  # Power simulation

  df = data.frame(power=c(unlist(lapply(a_til_pow,FUN=mean)),unlist(lapply(a_til_s_pow,FUN=mean)),unlist(lapply(a_til_p_pow,FUN=mean)),unlist(lapply(a_til_s_p_pow,FUN=mean)),unlist(lapply(a_hat_pow,FUN=mean))),estimate=c(rep(c(rep("z",5),rep("z*",5)),2),rep("MLE",5)),BC=c(rep("No",10),rep("Yes",10),rep("No",5)),SS=rep(n,5))

  power_fun = function(n){  
    X = cbind(rep(1,maxT),1:maxT)
    W = diag(exp(b0+b2[l]*(1:maxT))/(1+exp(b0+b2[l]*(1:maxT)))^2)
    se = sqrt(solve(t(X) %*% W %*% X)[2,2])
    return(sum(pnorm(abs(b2[l])/(se/sqrt(1:n))-qnorm(.975))*dbinom(1:n,n,t1)))
  }

  df2 = data.frame(power=unlist(lapply(n,power_fun)),SS=n)


  p=ggplot(df,aes(x=SS,y=power,color=estimate,linetype=BC)) + geom_line() + ylab("Power") + xlab("Sample Size") + labs(color="Estimate",linetype="Bias-Corrected") + geom_line(data=df2,mapping=aes(x=SS,y=power),color="Black",linetype="dashed")

  ggsave(paste("slope_power_fm_",k,".png",sep=""),p)

}
