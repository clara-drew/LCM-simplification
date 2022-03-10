library(dplyr)
library(parallel)
library(flexmix)
source("simulation_funs.R")

# Read in parameters for particular simulation parameters

source("parameters/s2b-beta01-parameters.R")

set.seed(122)

seeds <- list()
for(i in 1:length(maxT)){
  seeds[[i]] <- sample(1000:9999,size=its,replace=FALSE)
}  

# Create lists for results

theta_hat = list()
theta_hat_BC = list()
theta_hat_s = list()
theta_hat_s_BC = list()
theta_hat_MLE = list()

alpha_hat = list()
alpha_hat_BC = list()
alpha_hat_s = list()
alpha_hat_s_BC = list()
alpha_hat_MLE = list()

beta_hat_s= list()
beta_hat_s_BC = list()
beta_hat_MLE = list()

theta_hat_CP = list()
theta_hat_BC_CP = list()
theta_hat_s_CP = list()
theta_hat_s_BC_CP = list()
theta_hat_MLE_CP = list()

alpha_hat_CP = list()
alpha_hat_BC_CP = list()
alpha_hat_s_CP = list()
alpha_hat_s_BC_CP = list()
alpha_hat_MLE_CP = list()

beta_hat_s_CP= list()
beta_hat_s_BC_CP= list()
beta_hat_MLE_CP= list()

for(j in 3:4){ # Iterate simulation for differing numbers of tests
  
  alpha=alpha_list[[j]] # Set alpha(t) based on number of tests and parameter specifications
  
  # Initialize empty vectors/matrix for jth iteration
  
  theta_hat[[j]] = rep(NA,its)
  theta_hat_s[[j]] = rep(NA,its)
  alpha_hat[[j]] = matrix(ncol=maxT[j],nrow=its)
  alpha_hat_s[[j]] = matrix(ncol=maxT[j],nrow=its)
  theta_hat_BC[[j]] = rep(NA,its)
  theta_hat_s_BC[[j]] = rep(NA,its)
  alpha_hat_BC[[j]] = matrix(ncol=maxT[j],nrow=its)
  alpha_hat_s_BC[[j]] = matrix(ncol=maxT[j],nrow=its)
  theta_hat_MLE[[j]] = rep(NA,its)
  alpha_hat_MLE[[j]] = matrix(ncol=maxT[j],nrow=its)
  
  beta_hat_MLE[[j]] = rep(NA,its)
  beta_hat_s[[j]] = rep(NA,its)
  beta_hat_s_BC[[j]] = rep(NA,its)
  
  theta_hat_CP[[j]] = rep(NA,its)
  theta_hat_BC_CP[[j]] = rep(NA,its)
  theta_hat_s_CP[[j]] = rep(NA,its)
  theta_hat_s_BC_CP[[j]] = rep(NA,its)
  alpha_hat_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  alpha_hat_BC_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  alpha_hat_s_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  alpha_hat_s_BC_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  theta_hat_MLE_CP[[j]] = rep(NA,its)
  alpha_hat_MLE_CP[[j]] = matrix(ncol=maxT[j],nrow=its)
  
  beta_hat_MLE_CP[[j]] = rep(NA,its)
  beta_hat_s_CP[[j]] = rep(NA,its)
  beta_hat_s_BC_CP[[j]] = rep(NA,its)
  
  for(it in 1:its){ # iterate its times (specified in parameters file)
    
    set.seed(seeds[[j]][it])
    
    zeta=rbinom(n,1,theta)
    t=rep(0:(maxT[j]-1),n)
    ids=rep(1:n,each=maxT[j])
    y=rep(NA,maxT[j]*n)
    
    alphai <- alpha_fun(maxT[j])
    
    for(i in 1:n){
      for(l in 0:(maxT[j]-1)){
        if(zeta[i]==1){y[ids==i & t==l]=rbinom(1,1,alphai[i,l+1])}
        else{y[ids==i & t==l]=rbinom(1,1,beta)}
      }
    }
    
    # z, simple
    
    df = data.frame(t=t,ids=ids,y=y)
    z_samp = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>0,1,0)
    df_samp = left_join(data.frame(ids=1:n,z=z_samp),df)
    
    theta_samp=mean(z_samp)
    alpha_samp = aggregate(y~t,data=df_samp[df_samp$z==1,],FUN=mean)[,2]    
    theta_hat[[j]][it]=theta_samp
    alpha_hat[[j]][it,]=alpha_samp
    
    theta_hat_CP[[j]][it]=cover(theta,theta_samp,n)
    n1=sum(z_samp)
    alpha_hat_CP[[j]][it,]=cover(alpha,alpha_samp,n1)
    
    # z-star, simple
    
    z_samp1 = df[!duplicated(df$id),]$y
    z_samp2 = ifelse(aggregate(df$y,by=list(ids=df$ids),FUN=sum)[,2]>1,1,0)
    z_samp_s= ifelse(z_samp1+z_samp2>0,1,0)
    df_samp_s= left_join(data.frame(ids=1:n,z=z_samp_s),df)
    theta_samp_s = mean(z_samp_s)
    alpha_samp_s = aggregate(y~t,data=df_samp_s[df_samp_s$z==1,],FUN=mean)[,2]
    beta_samp_s = mean(df_samp_s[df_samp_s$z==0,]$y)
    theta_hat_s[[j]][it]=theta_samp_s
    alpha_hat_s[[j]][it,]=alpha_samp_s
    beta_hat_s[[j]][it]=beta_samp_s
    
    theta_hat_s_CP[[j]][it]=cover(theta,theta_samp_s,n)
    n1=sum(z_samp_s)
    n0=dim(df_samp_s[df_samp_s$z==0,])[1]
    alpha_hat_s_CP[[j]][it,]=cover(alpha,alpha_samp_s,n1)
    beta_hat_s_CP[[j]][it]=cover(beta,beta_samp_s,n0)
    
    # z, bias corrected
    
    gp = rep(NA,n)
    
    for(i in 1:n){
      tps = df_samp[df_samp$ids==i,]$t+1
      gp[i]=prod(1-alpha_samp[tps],na.rm=T)  
    }
    
    thetaHatBC=(theta_samp)/(mean(1-gp))
    theta_hat_BC[[j]][it]=thetaHatBC
    alpha_hat_BC[[j]][it,]=(alpha_samp*mean(1-gp))
    
    boot_ci = bootCI(its=1000,cores=ncores)
    theta_hat_BC_CP[[j]][it]=coverp(theta,boot_ci,1)
    alpha_hat_BC_CP[[j]][it,]=coverp(alpha,boot_ci,2:(1+maxT[j]))
    
    # z-star, bias corrected
    
    pz10 = rep(NA,n)
    pz11 = rep(NA,n)
    pz111 = matrix(nrow=n,ncol=maxT[j])
    pz110 = rep(NA,n)
    
    for(i in 1:n){
      tps = df_samp_s[df_samp_s$ids==i,]$t+1
      Ti = length(tps)
      pz10[i] = 1-sum(dbinom(0:1,Ti-1,beta_samp_s)*(1-beta_samp_s))
      pz11[i] = 1-prod(1-alpha_samp_s[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) alpha_samp_s[x]*prod(1-alpha_samp_s[-x]))))      
      pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-alpha_samp_s[tps[-x]]))))
      pz110[i] = 1-(1-beta_samp_s)^(Ti-1)
    }
    
    thetaHatBC=(theta_samp_s-mean(pz10))/mean(pz11-pz10)
    theta_hat_s_BC[[j]][it]=thetaHatBC
    alphaHatBC=(alpha_samp_s*(mean(pz11)*thetaHatBC+mean(pz10)*(1-thetaHatBC))-
              mean(pz110)*beta_samp_s*thetaHatBC)/(colMeans(pz111,na.rm=T)*thetaHatBC)
    alpha_hat_s_BC[[j]][it,]=alphaHatBC
    beta_hat_s_BC[[j]][it]=(beta_samp_s*(mean(1-pz11)*thetaHatBC+mean(1-pz10)*(1-thetaHatBC))-
                           mean(colMeans(1-pz111,na.rm=T)*alphaHatBC*(1-thetaHatBC)))/mean((1-pz110)*thetaHatBC)
    
    boot_ci= bootCI_s(its=1000,cores=ncores)
    
    theta_hat_s_BC_CP[[j]][it]=coverp(theta,boot_ci,1)
    alpha_hat_s_BC_CP[[j]][it,]=coverp(alpha,boot_ci,2:(1+maxT[j]))
    beta_hat_s_BC_CP[[j]][it]=coverp(beta,boot_ci,2+maxT[j])
    
    # EM algorithm MLE using flexmix function
    
    df <- data.frame(y,t,ids)
    
    df_short <- df %>% group_by(ids) %>% filter(row_number() == 1)
    df$Ez <- NA
    
    k=1
    for(i in 1:dim(df_short)[1]){
      times<-as.numeric(table(df$ids)[i])
      df$Ez[k:(k+times-1)] <- rep(df_short$y[i]+1,times)
      k<-k+times
    }
    
    cluster <- df$Ez
    
    if(mean(df$Ez)==1){
      cluster <- NULL
    }
    
    m1 <- flexmix(cbind(y,1-y) ~ 1 | ids, data=df, k=2, cluster=cluster,
                  model = FLXMRglmfix(nested = list(k=c(1,1), formula = c(~0,~t)),family="binomial"))
    
    refit_m1 <- refit(m1)
    
    coefs <- refit_m1@coef
    Cov <- refit_m1@vcov
    SEs <- sqrt(diag(Cov))
    
    if(length(coefs)!=1){ # coefs will be length 1 if converged to single shedding group

      theta_hat_MLE[[j]][it] = invlogit(coefs[4])
      alpha_hat_MLE[[j]][it,] = invlogit(coefs[3]+coefs[1]*(0:(maxT[j]-1)))
      beta_hat_MLE[[j]][it] = invlogit(coefs[2])
    
      ME1 = qnorm(.975)*SEs[4]
      theta_hat_MLE_CP[[j]][it] = as.numeric(logit(theta)>coefs[4]-ME1 & logit(theta)<coefs[4]+ME1)
      for(t in 1:maxT[j]){
        dqdb = c((t-1),0,1,0)
        ME2 = qnorm(.975)*sqrt(dqdb %*% Cov %*% dqdb)
        alpha_hat_MLE_CP[[j]][it,t] = as.numeric(logit(alpha[t])>coefs[3]+coefs[1]*(t-1)-ME2 & logit(alpha[t])<coefs[3]+coefs[1]*(t-1)+ME2)
      }
      ME3 = qnorm(.975)*SEs[2]
      beta_hat_MLE_CP[[j]][it] = as.numeric(logit(beta[j])>coefs[2]-ME3 & logit(beta[j])<coefs[2]+ME3)
    }
    
    # Save results in external files (folder name specified in parameters file)
    
    df_theta = data.frame(theta_hat=theta_hat[[j]],theta_hat_BC=theta_hat_BC[[j]],theta_hat_s=theta_hat_s[[j]],theta_hat_s_BC=theta_hat_s_BC[[j]],theta_hat_MLE=theta_hat_MLE[[j]])
    df_alpha = cbind(alpha_hat[[j]],alpha_hat_BC[[j]],alpha_hat_s[[j]],alpha_hat_s_BC[[j]],alpha_hat_MLE[[j]])
    colnames(df_alpha)=c(paste("alpha_hat_",1:maxT[j],sep=""),paste("alpha_hat_BC_",1:maxT[j],sep=""),paste("alpha_hat_s_",1:maxT[j],sep=""),paste("alpha_hat_s_BC_",1:maxT[j],sep=""),paste("alpha_hat_MLE_",1:maxT[j],sep=""))
    write.csv(df_theta,paste(folder,"/df_","theta_",maxT[j],".csv",sep=""))
    write.csv(df_alpha,paste(folder,"/df_","alpha_",maxT[j],".csv",sep=""))
    df_theta_cp= data.frame(theta_hat=theta_hat_CP[[j]],theta_hat_BC=theta_hat_BC_CP[[j]],theta_hat_s=theta_hat_s_CP[[j]],theta_hat_s_BC=theta_hat_s_BC_CP[[j]],theta_hat_MLE=theta_hat_MLE_CP[[j]])
    df_alpha_cp= cbind(alpha_hat_CP[[j]],alpha_hat_BC_CP[[j]],alpha_hat_s_CP[[j]],alpha_hat_s_BC_CP[[j]],alpha_hat_MLE_CP[[j]])
    colnames(df_alpha_cp)=c(paste("alpha_hat_",1:maxT[j],sep=""),paste("alpha_hat_BC_",1:maxT[j],sep=""),paste("alpha_hat_s_",1:maxT[j],sep=""),paste("alpha_hat_s_BC_",1:maxT[j],sep=""),paste("alpha_hat_MLE_",1:maxT[j],sep=""))
    write.csv(df_theta_cp,paste(folder,"/df_","theta_cp_",maxT[j],".csv",sep=""))
    write.csv(df_alpha_cp,paste(folder,"/df_","alpha_cp_",maxT[j],".csv",sep=""))
    

  }
  
}