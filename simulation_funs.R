# Indicator function for c95% confidence interval coverage
# for sample proportion

cover <- function(p,phat,n){
  me = sqrt(phat*(1-phat)/n)
  return(ifelse(p>phat-qnorm(.975)*me & p<phat+qnorm(.975)*me,1,0))
}

# Indicator function for bootstrapped 95% confidence interval coverage
# for bias corrected sample proportion 

coverp <- function(p,boot_ci,i){
  cover=c()
  for(k in 1:length(i)){
    cover=c(cover,ifelse(p[k]>boot_ci[[i[k]]][1] & p[k]<boot_ci[[i[k]]][2],1,0))
  }
  return(cover)
}

# Function for single bootstrap iteration for bias correction
# of z indicator estimates

iter <- function(x){
  ids=sample(1:n,n,replace=T)
  dfB=data.frame(t=NULL,ids=NULL,y=NULL)
  for(i in 1:n){
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df[df$ids==ids[i],])
    ni = dim(df[df$ids==ids[i],])[1]
    dfB$ids[(nl+1):(nl+ni)]=i
  }
  
  z_samp = ifelse(aggregate(dfB$y,by=list(ids=dfB$ids),FUN=sum)[,2]>0,1,0)
  
  if(mean(z_samp)==0){
    return(list(theta_hat_BC=NA,alpha_hat_BC=rep(NA,max(dfB$t)+1)))
  }
  
  else{
  
    ids = 1:n
    df_boot = left_join(data.frame(ids=ids,z=z_samp),dfB)
  
    gp = rep(NA,length(ids))
  
    theta_hat = mean(z_samp)

    alpha_hat = aggregate(y~t,data=df_boot[df_boot$z==1,],FUN=mean)[,2]
  
    for(i in 1:n){
      tps = df_boot[df_boot$ids==i,]$t+1
      ntps = length(tps)
      gp[i]=prod(1-alpha_hat[tps],na.rm=T)  
    }
  
    theta_hat_BC=(theta_hat)/(mean(1-gp))
    alpha_hat_BC=(alpha_hat*sum((1-gp)*theta_hat_BC))/(theta_hat_BC*n)
  
    return(list(theta_hat_BC=theta_hat_BC,alpha_hat_BC=alpha_hat_BC))
  }  
}

# Function for percentile bootstrapped 95% confidence interval for
# bias corrected z indicator estimates

bootCI <- function(its=1000,cores=3){
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

# Function for single bootstrap iteration for bias correction
# of z-star indicator estimates

iter_s <- function(x){
  
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
  
  if(mean(z_samp)==0){
    return(list(theta_hat_s_BC=NA,alpha_hat_s_BC=rep(NA,max(dfB$t)+1),beta_hat_s_BC=NA))
  }
  
  else{
  
    df_boot = left_join(data.frame(ids=1:n,z=z_samp),dfB)
    theta_hat_s=mean(z_samp)
    alpha_hat_s = aggregate(df_boot[df_boot$z==1,]$y,by=list(ids=df_boot[df_boot$z==1,]$t),FUN=mean)[,2]
    beta_hat_s = mean(df_boot[df_boot$z==0,]$y)
  
    pz10 = rep(NA,n)
    pz11 = rep(NA,n)
    pz111 = matrix(nrow=n,ncol=maxT[j])
    pz110 = rep(NA,n)
  
    for(i in 1:n){
      tps = df_boot[df_boot$ids==i,]$t+1
      Ti = length(tps)
      pz10[i] = (1-sum(dbinom(0:1,Ti-1,beta_hat_s)*(1-beta_hat_s)))
      pz11[i] = (1-prod(1-alpha_hat_s[tps])-sum(unlist(lapply(2:Ti,FUN=function(x) alpha_hat_s[x]*prod(1-alpha_hat_s[-x])))))       
      pz111[i,tps] = c(1,unlist(lapply(tps[-1],FUN = function(x) 1-prod(1-alpha_hat_s[tps[-x]]))))
      pz110[i] = 1-(1-beta_hat_s)^(Ti-1)
    }
  
    theta_hat_s_BC=(theta_hat_s-mean(pz10))/mean(pz11-pz10)
    alpha_hat_s_BC=(alpha_hat_s*(mean(pz11)*theta_hat_s_BC+mean(pz10)*(1-theta_hat_s_BC))-
            mean(pz110)*beta_hat_s*theta_hat_s_BC)/(colMeans(pz111,na.rm=T)*theta_hat_s_BC)
    beta_hat_s_BC=(beta_hat_s*(mean(1-pz11)*theta_hat_s_BC+mean(1-pz10)*(1-theta_hat_s_BC))-
            mean(colMeans(1-pz111,na.rm=T)*alpha_hat_s_BC*(1-theta_hat_s_BC)))/mean((1-pz110)*theta_hat_s_BC)
  
  
    return(list(theta_hat_s_BC=theta_hat_s_BC,alpha_hat_s_BC=alpha_hat_s_BC,beta_hat_s_BC=beta_hat_s_BC))

    }
}

# Function for percentile bootstrapped 95% confidence interval for
# bias corrected z-star indicator estimates

bootCI_s <- function(its=1000,cores=3){
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

# Logit and inverse logit functions

invlogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# Logistic alpha(ti) using specified parameters in file (scenarios 1 and 2)

alpha_logistic <- function(maxT){
  rand <- rnorm(n,0,sigma)
  alphai <- matrix(nrow=n,ncol=maxT)
  for(i in 1:n){
    alphai[i,] <- invlogit(a0+a2[j]*(0:(maxT-1))+rand[i])
  }
  return(alphai)
}

# Logistic mixture alpha(ti) using specified parameters in file (scenario 2b)

alpha_logistic_unident <- function(maxT){
  rand <- rnorm(n,0,sigma)
  alphai <- matrix(nrow=n,ncol=maxT)
  for(i in 1:n){
    vit <- a0+a2[j]*(0:(maxT-1))+rand[i]
    alphai[i,] <- invlogit(omega0+omega1*vit)
  }
  return(alphai)
}

# Helper function for polynomial alpha

poly <- function(x,cov){
  p <- c(x^4,x^3,x^2,x,1)
  return(t(cov) %*% p)
}

# Polynomial alpha(t) using specified parameters in file (scenario 3)

alpha_polynomial <- function(maxT){
  min <- 19/3
  last <- 1
  y1 <- 0.875
  y2 < -0.3
  y3 <- 0.125
  A = rbind(c(0,0,0,1,0),c(4*min^3,3*min^2,2*min,1,0),c(0,0,0,0,1),
            c(min^4,min^3,min^2,min,1),c(last^4,last^3,last^2,last,1))
  b = c(0,0,y1,y2,y3)
  cov <- solve(A) %*% b
  tps <- ((1:maxT)-1)/(maxT-1)
  one_row <- unlist(lapply(tps,poly,cov=cov))
  return(matrix(nrow=n,ncol=length(one_row),data=one_row,by=2))
}


