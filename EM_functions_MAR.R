# Load packages

library(ggplot2)
library(parallel)
library(plyr)
library(gridExtra)

# Functions

invlogit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

# log-likelihood function to find beta0 and beta2

loglik=function(x,z,y,t){
  b0=x[1]
  b2=x[2]
  return(sum(z*y*b0+z*y*b2*t-z*log(1+exp(b0+b2*t)),na.rm=T))
}

# Expected Value Functions

Ez=function(b0,b1,b2,b3,y,t){ # Takes y_t vector for person i, returns single value
  obs=which(!is.na(y))
  yi=y[obs]
  ti=t[obs]
  n=length(yi)
  (exp(b1)*prod(invlogit(b0+b2*ti)^yi)*prod((1-invlogit(b0+b2*ti))^(1-yi)))/
    (exp(b1)*prod(invlogit(b0+b2*ti)^yi)*prod((1-invlogit(b0+b2*ti))^(1-yi))+
       invlogit(b3)^sum(yi)*(1-invlogit(b3))^(n-sum(yi)))
}


# Generate Random Data

gen_dat = function(beta0,beta1,beta2,beta3,t,n){
  # Person id and time, t
  df=data.frame(id=rep(1:n,each=t),t=rep(1:t,n))
  # Probabilities of shedding at each time point given person is shedder
  theta2=invlogit(beta0+beta2*(1:t))
  # Generate shedder/nonshedder (zeta_i) data from each person
  df$z=rep(rbinom(n,1,invlogit(beta1)),each=t)
  # Generate test results given zeta
  df$y=NA
  for(i in 1:t){
    len=length(df[df$t==i & df$z==1,]$y)
    df[df$t==i & df$z==1,]$y=rbinom(len,1,theta2[i])
  }
  len=length(df[df$z==0,]$z)
  df[df$z==0,]$y=rbinom(len,1,invlogit(beta3))
  # Generate data for whether test result is missing at time, t, given result of previous test
  df$r=NA

  first = trunc(rtruncnorm(n,a=1,b=51,mean=20,sd=8))

  for(i in 1:n){
    if(first[i]>1){
      df$y[(1+t*(i-1)):(t*(i-1)+first[i]-1)]=NA
    } 
    phis=ifelse(df$y[(t*(i-1)+first[i]):(t*i-1)]==1,phi1,phi2)
    if(first[i]<=(t-1)){
      df$r[(t*(i-1)+first[i]+1):(t*i)]=rbinom(50-first[i],1,phis)
    }  
  }  
   

  df[df$r==0,]$y=NA
  
  
  return(df)
}  

# Impute expected value

impute=function(df,b0,b1,b2,b3){ # Returns data set with current expected values imputed

  # Create variable for expected value of zeta using function, Ez
  if(min(df$t)==0){
    t=max(df$t)+1
  }
  else{
    t=max(df$t)
  }  
  df$ez=NA
  n=length(unique(df$id))
  for(i in 0:(n-1)){
    df$ez[(t*i+1):(t*i+t)]=Ez(b0,b1,b2,b3,y=df$y[(t*i+1):(t*i+t)],t=df$t[(t*i+1):(t*i+t)])
  }
  
  return(df)
  
}  

# Find Estimates

next_b1 = function(ez){
  logit(mean(ez))
}

next_b3 = function(ez,y){
  obs = which(!is.na(y))
  y=y[obs]
  ez=ez[obs]
  logit((sum(y-ez*y))/(sum(1-ez)))
}

# Use optimize and log-likelihood function to find next estimate
# for beta0 and beta2. Use previous estimates as first guess.

next_b02 = function(b0,b2,z,y,t){
  x=c(b0,b2)
  obs = which(!is.na(y))
  y=y[obs]
  z=z[obs]
  t=t[obs]
  out=optim(par=x,fn=loglik,z=z,y=y,t=t,control=list(fnscale=-1))
  return(list(b0=out$par[1],b2=out$par[2]))
}

# Function for one iteration of simulation given beta values, 
# initial guesses, time points and number of participants

iter=function(inits, betas, t, n){
  
  b0est=NA
  b1est=NA
  b2est=NA
  b3est=NA
  
  beta0=betas[1]
  beta1=betas[2]
  beta2=betas[3]
  beta3=betas[4]
  
  b0_init=inits[1]
  b1_init=inits[2]
  b2_init=inits[3]
  b3_init=inits[4]
  
  df=gen_dat(beta0,beta1,beta2,beta3,t,n)

  b0est[1]=b0_init
  b1est[1]=b1_init
  b2est[1]=b2_init
  b3est[1]=b3_init

  k=1

  df=impute(df,b0_init,b1_init,b2_init,b3_init)

  k=2

  nextb02 = next_b02(b0est[1],b2est[1],df$ez,df$y,df$t)

  b0est[2]=nextb02$b0
  b1est[2]=next_b1(df$ez)
  b2est[2]=nextb02$b2
  b3est[2]=next_b3(df$ez,df$y)
  
  # Continue finding next estimate until estimates are within 0.001
  
  while(abs(b0est[k]-b0est[k-1])>0.001 | abs(b1est[k]-b1est[k-1])>0.001 | 
    abs(b2est[k]-b2est[k-1])>0.001 | abs(b3est[k]-b3est[k-1])>0.001){
    
    df=impute(df,b0est[k],b1est[k],b2est[k],b3est[k])
    
    k=k+1

    nextb02 = next_b02(b0est[k-1],b2est[k-1],df$ez,df$y,df$t)
    
    b0est[k]=nextb02$b0
    b1est[k]=next_b1(df$ez)
    b2est[k]=nextb02$b2
    b3est[k]=next_b3(df$ez,df$y)
    
  }
  
  return(list(b0est=b0est[k],b1est=b1est[k],b2est=b2est[k],b3est=b3est[k],k=k))
}

getval = function(res,its,k){
  return(res[[its]][k])
}

# Function to run iterations in parallel and return estimates for each parameter

sim=function(niters,inits,betas,t,n){
  seed=seq(from=66,to=66+niters-1)
  res=mclapply(1:niters,function(it){
    set.seed(seed[it])
    result=try(iter(inits,betas,t,n),silent = T)
    if(!("try-error" %in% class(result))){
      return(c(result$b0est,result$b1est,result$b2est,result$b3est,max(result$k),seed[it]))
    }
    else{
      return(rep(NA,6))
    }
  },mc.cores=12)
  b0est=unlist(lapply(1:niters,getval,res=res,k=1))
  b1est=unlist(lapply(1:niters,getval,res=res,k=2))
  b2est=unlist(lapply(1:niters,getval,res=res,k=3))
  b3est=unlist(lapply(1:niters,getval,res=res,k=4))
  its=unlist(lapply(1:niters,getval,res=res,k=5))
  seed=unlist(lapply(1:niters,getval,res=res,k=6))
  df=data.frame(b0est,b1est,b2est,b3est,its,seed)
  return(df)
}

# Get phi estimates

get_phi = function(df){
  ids = unique(df$id)

  t = max(df$t)
  n = length(unique(df$id))


  min_t=rep(NA,length(ids))
  max_t=rep(NA,length(ids))

  for(i in 1:n){
    min_t[i] = min(df[!is.na(df$y) & df$id==ids[i],]$t)
    max_t[i] = max(df[!is.na(df$y) & df$id==ids[i],]$t)
  }


  y1=0
  y0=0
  r1=0
  r0=0

  for(i in 1:length(min_t)){
    k=min_t[i]
    y=df[df$id==ids[i],]$y
    while(k<max_t[i]){
      if(y[k]==1){
        k=k+1
        while(is.na(y[k])){
          y1=y1+1
          k=k+1
          if(k>max_t[i]-1) break
        }
        if(k>max_t[i]-1) break
        y1=y1+1
        r1=r1+1
      }
      else{
        k=k+1
        while(is.na(y[k])){
          y0=y0+1
          k=k+1
          if(k>max_t[i]-1) break
        }
        if(k>max_t[i]-1) break
        y0=y0+1
        r0=r0+1
      }
    }    
  }


  p1est=r1/y1
  p0est=r0/y0

  return(list(p1est=p1est,p0est=p0est))
}

get_init = function(df){
  ids=unique(df$id)
  z=NULL
  for(i in 1:length(ids)){
    df_id=df[df$id==ids[i],]
    obs=which(!is.na(df_id$y))
    if(length(obs)==0){
      z=c(z,rep(NA,max(df$t)))
    }
    else{  
      first_t=df_id[min(which(!is.na(df_id$y))),]$t
      if(df[df$t==first_t & df$id==ids[i],]$y==1){
        z=c(z,rep(1,length(df[df$id==ids[i],]$y)))
      }  
      else{
        z=c(z,rep(0,length(df[df$id==ids[i],]$y)))
      }
    }   
  }
  m1=tryCatch(glm(y~t,data=df[z==1,],family="binomial"),warning=function(x) return(NA))
  m2=tryCatch(glm(y~1,data=df[z==0,],family="binomial"),warning=function(x) return(NA))
  if(is.na(m1)){
    b0_init=3
    b2_init=-1
  }
  else{
    b0_init=coefficients(m1)[1]
    b2_init=coefficients(m1)[2]
  }
  if(is.na(m2)){
    b3_init=0.01
  }
  else{
    b3_init=coefficients(m2)[1]
  }
  b1_init=logit(mean(z,na.rm=T))
  return(list(b0_init=b0_init,b1_init=b1_init,b2_init=b2_init,b3_init=b3_init))
}

get_EM_est = function(df,tol=0.001){

  if(min(df$t)==0){
    t=max(df$t)+1
  }
  else{
    t=max(df$t)
  }
  n=length(unique(df$id))

  b0est=NA
  b1est=NA
  b2est=NA
  b3est=NA

  inits=get_init(df)

  b0est[1]=inits$b0_init
  b1est[1]=inits$b1_init
  b2est[1]=inits$b2_init
  b3est[1]=inits$b3_init

  k=1

  df=impute(df,b0est[1],b1est[1],b2est[1],b3est[1])

  k=2

  nextb02=next_b02(b0est[1],b2est[1],df$ez,df$y,df$t)

  b0est[2]=nextb02$b0
  b1est[2]=next_b1(df$ez)
  b2est[2]=nextb02$b2
  b3est[2]=next_b3(df$ez,df$y)

  while((abs(b0est[k]-b0est[k-1])>tol | abs(b1est[k]-b1est[k-1])>tol | abs(b2est[k]-b2est[k-1])>tol | abs(b3est[k]-b3est[k-1])>tol) & k<100){
  
    df=impute(df,b0est[k],b1est[k],b2est[k],b3est[k])
  
    k=k+1

    nextb02 = next_b02(b0est[k-1],b2est[k-1],df$ez,df$y,df$t)

    b0est[k]=nextb02$b0
    b1est[k]=next_b1(df$ez)
    b2est[k]=nextb02$b2
    b3est[k]=next_b3(df$ez,df$y)

    if(b3est[k]==-Inf) break

  }

  b0=b0est[k]
  b1=b1est[k]
  b2=b2est[k]
  b3=b3est[k]

  
  return(list(b0=b0,b1=b1,b2=b2,b3=b3,z=df$ez[seq(from=1,to=n*t,by=t)]))
}

calc_Cov_no3 = function(b0,b1,b2,z,t,y,id,n=265){

  ts=aggregate(t,list(id=id),FUN=length)[,2]

  varz=z*(1-z)

  # IX

  IX00 = sum(z*aggregate(invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IX01 = 0
  IX02 = sum(z*aggregate(t*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IX11 = n*exp(b1)/(1+exp(b1))^2
  IX12 = 0
  IX22 = sum(z*aggregate(t^2*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])

  IX = cbind(c(IX00,IX01,IX02),c(IX01,IX11,IX12),c(IX02,IX12,IX22))

  # IX|Y

  IXY00 = sum(varz*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2]^2)
  IXY01 = sum(varz*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2])
  IXY02 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2]*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2])
  IXY11 = sum(varz)
  IXY12 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IXY22 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2]^2)

  IXY = cbind(c(IXY00,IXY01,IXY02),c(IXY01,IXY11,IXY12),c(IXY02,IXY12,IXY22))

  FI = IX - IXY
  Cov = solve(FI,tol=1e-20)
  return(list(Cov=Cov,FI=FI,IX=IX,IXY=IXY))

}



calc_Cov = function(b0,b1,b2,b3,z,t,y,id,n=265){
 
  ts=aggregate(t,list(id=id),FUN=length)[,2]

  varz=z*(1-z)

  # IX

  IX00 = sum(z*aggregate(invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IX01 = 0
  IX02 = sum(z*aggregate(t*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IX03 = 0
  IX11 = n*exp(b1)/(1+exp(b1))^2
  IX12 = 0
  IX13 = 0
  IX22 = sum(z*aggregate(t^2*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IX23 = 0
  IX33 = sum((1-z)*(ts*invlogit(b3)*(1-invlogit(b3))))

  IX = cbind(c(IX00,IX01,IX02,IX03),c(IX01,IX11,IX12,IX13),c(IX02,IX12,IX22,IX23),c(IX03,IX13,IX23,IX33))

  # IX|Y

  IXY00 = sum(varz*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2]^2)
  IXY01 = sum(varz*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2])
  IXY02 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2]*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2])
  IXY03 = -sum(varz*aggregate(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t),by=list(id=id),FUN=sum)[,2]*aggregate(y*(1-invlogit(b3))-(1-y)*invlogit(b3),by=list(id=id),FUN=sum)[,2])
  IXY11 = sum(varz)
  IXY12 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2])
  IXY13 = -sum(varz*aggregate(y*(1-invlogit(b3))-(1-y)*invlogit(b3),by=list(id=id),FUN=sum)[,2])
  IXY22 = sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2]^2)
  IXY23 = -sum(varz*aggregate(t*(y*(1-invlogit(b0+b2*t))-(1-y)*invlogit(b0+b2*t)),by=list(id=id),FUN=sum)[,2]*aggregate(y*(1-invlogit(b3))-(1-y)*invlogit(b3),by=list(id=id),FUN=sum)[,2])
  IXY33 = sum(varz*aggregate(y*(1-invlogit(b3))-(1-y)*invlogit(b3),by=list(id=id),FUN=sum)[,2]^2)
      

  IXY = cbind(c(IXY00,IXY01,IXY02,IXY03),c(IXY01,IXY11,IXY12,IXY13),c(IXY02,IXY12,IXY22,IXY23),c(IXY03,IXY13,IXY23,IXY33))

  FI = IX - IXY
  Cov = solve(FI,tol=1e-20)
  return(list(Cov=Cov,FI=FI,IX=IX,IXY=IXY))
}

dqdb = function(b0,b2,t) return(c(invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),0,t*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),0))

dpdb = function(b1){
  return(c(0,invlogit(b1)*(1-invlogit(b1)),0,0))
}

dgdb = function(b0,b1,b2,t){
  J = rbind(c(0,invlogit(b1)*(1-invlogit(b1)),0,0),dqdb(b0,b2,t))
  return(J)
}

get_set1=function(Cov,b1){
  return(sqrt(dpdb(b1) %*% Cov %*% dpdb(b1)))
}

get_set2_t = function(t,Cov,b0,b2){
  dtdb0 = invlogit(b0+b2*t)*(1-invlogit(b0+b2*t))
  dtdb = c(dtdb0,t*dtdb0)
  return(sqrt(dtdb %*% Cov %*% dtdb))
}

get_set1t2_t = function(t,Cov,b0,b1,b2){
  t1=invlogit(b1)
  t2=invlogit(b0+b2*t)
  C = dgdb(b0,b1,b2,t) %*% Cov %*% t(dgdb(b0,b1,b2,t))
  Var = 4*t1*t2*C[1,2]+2*C[1,2]^2+(C[1,1]+t1^2)*(C[2,2]+t2^2)-(C[1,2]+t1*t2)^2
  return(sqrt(Var))
}

get_set1t2_tno3 = function(t,Cov,b0,b1,b2){
  t1=invlogit(b1)
  t2=invlogit(b0+b2*t)
  C = dgdb(b0,b1,b2,t)[,1:3] %*% Cov %*% t(dgdb(b0,b1,b2,t)[,1:3])
  Var = 4*t1*t2*C[1,2]+2*C[1,2]^2+(C[1,1]+t1^2)*(C[2,2]+t2^2)-(C[1,2]+t1*t2)^2
  return(sqrt(Var))
}

se_t1t2_simple = function(t,b0,b2,t1,cov,t2=NA){
  if(is.na(t2)){
     t2=invlogit(b0+b2*t)   
     # Find derivative of theta_2(t)
     dqdb0 = c(invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)))
     dqdb = c(dqdb0,t*dqdb0)
     # Delta method for variance of theta_2(t)
     var_t2 = dqdb%*%cov%*%dqdb
     # Variance of t1
  }
  else{var_t2=t2*(1-t2)/n}
  var_t1 = t1*(1-t1)/n
  Var = t1^2*var_t2+t2^2*var_t1+var_t1*var_t2
  return(sqrt(Var))

}

get_simp_est = function(df,logit=TRUE){
  n = length(unique(df$id))
  maxT = max(df$t)
  t1 = length(unique(df[df$y==1,]$id))/n
  zest = rep(ifelse(aggregate(df$y, by=list(id=df$id),FUN=sum,na.rm=T)$x>0,1,0),each=maxT)
  if(logit==TRUE){
    m = glm(y~t,family="binomial",data=df[zest==1,])
    cov = vcov(m)
    b0 = coef(m)[1]
    b2 = coef(m)[2]
    return(list(b0=b0,b2=b2,t1=t1,cov=cov))
  }
  else{
    t2=aggregate(df[zest==1,]$y, by=list(t=df[zest==1,]$t),FUN=mean,na.rm=T)[,2]
    return(list(t1=t1,t2=t2))
  }
}



bias_func = function(df,b0,t1,b2,t3,tps){
  
  pids = unique(df$id)
  N = length(pids)
  obs = !is.na(df$y)
  id = df$id[obs]
  t = df$t[obs]
  t2=invlogit(b0+b2*tps)
  gp = rep(NA,N)
  hp = rep(NA,N)

  for(i in 1:N){
    ti = t[id==pids[i]]
    gp[i]=prod(1-invlogit(b0+b2*ti))
    hp[i]=(1-t3)^(length(ti))
  }

  bias = mean(-gp)*t1*t2 + (1-t1)*t3*mean(1-hp)

  return(bias)
}

get_bias_corr_est2 = function(df,b0_til,t1_til,b2_til,t3,tps,t2_til=NULL){
  pids = unique(df$id)
  N = length(pids)
  obs = !is.na(df$y)
  id = df$id[obs]
  t = df$t[obs]
  gp = rep(NA,N)
  hp = rep(NA,N)
  pz = rep(NA,N)
  if(is.null(t2_til)){
    t2_til=invlogit(b0_til+b2_til*tps)
    for(i in 1:N){
      ti = t[id==pids[i]]
      gp[i]=prod(1-invlogit(b0_til+b2_til*ti))
      hp[i]=(1-t3)^(length(ti))
    }
  }

  else{
    for(i in 1:N){
      ti = t[id==pids[i]]
      gp[i]=prod(1-t2_til[ti])
      hp[i]=(1-t3)^(length(ti))
    }
  }

  t1 = (t1_til - mean(1-hp)) / (mean(-gp)+mean(hp))
  pzi = (1-gp)*t1 + (1-hp)*(1-t1)
  pzetai = ((1-gp)*t1) / ((1-gp)*t1+(1-hp)*(1-t1))
  t2 = ( t2_til*sum(pzi) - t3*sum((1-pzetai)*pzi) ) / sum(pzetai*pzi)

  return(list(t1=t1,t2=t2,t1t2=t1*t2))
}

#get_bias_corr_est = function(df,b0,t1,b2,t3,tps){
#  t2 = invlogit(b0+b2*tps)
#  return(t1*t2 - bias_func(df,b0,t1,b2,t3,tps))
#}

bootstrap_fun = function(its,df,t3,cores,logit=TRUE){
  ids = unique(df$id)
  N = length(ids)
  t = max(df$t)
  tps = unique(df$t)
  #t1t2 = matrix(ncol=t,nrow=its)
  t1t2_list= mclapply(1:its,function(x){
    ids_star = sample(ids,N,replace=T)
    df_star = df
    for(j in 1:N){
      df_star[(1+(j-1)*t):(j*t),] = df[df$id==ids_star[j],]
    }
    df_star$id = df$id
    ests_star = get_simp_est(df_star,logit)
    return(get_bias_corr_est2(df=df_star,b0_til=ests_star$b0,t1_til=ests_star$t1,b2_til=ests_star$b2,t3=t3,tps=tps,t2_til=ests_star$t2))
    },mc.cores=cores)
 
  #for(i in 1:its){
  #  ids_star = sample(ids,N,replace=T)
  #  df_star = df
  #  for(j in 1:N){
  #    df_star[(1+(j-1)*t):(j*t),] = df[df$id==ids_star[j],]
  #  }
  #  df_star$id = df$id
  #  ests_star = get_simp_est(df_star)
  #  t1t2[i,]=get_bias_corr_est(ests_star$b0,ests_star$t1,ests_star$b2,t3,tps)
  #}
  t1 = unlist(t1t2_list)[c(TRUE,rep(FALSE,100))]
  t2 = matrix(unlist(t1t2_list)[c(FALSE,rep(TRUE,50),rep(FALSE,50))],nrow=its,byrow=T)
  t1t2 = matrix(unlist(t1t2_list)[c(FALSE,rep(FALSE,50),rep(TRUE,50))],nrow=its,byrow=T)
  return(list(t1=t1,t2=t2,t1t2=t1t2))
}

get_CI_bias_corr= function(df,t3,its,cores,logit=TRUE){
  bootstraps = bootstrap_fun(its,df,t3,cores,logit)
  t1 = bootstraps$t1
  t2 = bootstraps$t2
  t1t2 = bootstraps$t1t2
  t2_list= as.list(as.data.frame(t2))
  t1t2_list= as.list(as.data.frame(t1t2))

  lower_t1 = quantile(t1,probs=0.025)
  upper_t1 = quantile(t1,probs=0.975)
  CIs_t2 = lapply(t2_list,quantile,probs=c(0.025,0.975))
  lower_t2 = unlist(CIs_t2)[c(TRUE,FALSE)]
  upper_t2 = unlist(CIs_t2)[c(FALSE,TRUE)]
  CIs_t1t2= lapply(t1t2_list,quantile,probs=c(0.025,0.975))
  lower_t1t2= unlist(CIs_t1t2)[c(TRUE,FALSE)]
  upper_t1t2= unlist(CIs_t1t2)[c(FALSE,TRUE)]
  return(list(lower_t1=lower_t1,upper_t1=upper_t1,lower_t2=lower_t2,upper_t2=upper_t2,lower_t1t2=lower_t1t2,upper_t1t2=upper_t1t2))
}

bootstrap_fun2 = function(its,df,t3,cores,logit=TRUE){
  ids = unique(df$id)
  N = length(ids)
  t = max(df$t)
  tps = unique(df$t)
  #t1t2 = matrix(ncol=t,nrow=its)
  t1t2_list= mclapply(1:its,function(x){
    ids_star = sample(ids,N,replace=T)
    df_star = df
    for(j in 1:N){
      df_star[(1+(j-1)*t):(j*t),] = df[df$id==ids_star[j],]
    }
    df_star$id = df$id
    ests_star = get_simp_est(df_star,logit)
    ests_bc = get_bias_corr_est2(df=df_star,b0_til=ests_star$b0,t1_til=ests_star$t1,b2_til=ests_star$b2,t3=t3,tps=tps,t2_til=ests_star$t2)
    n1 = round(N*ests_bc$t1)
    counts = round(n1*ests_bc$t2)
    y1 = unlist(lapply(counts,FUN=function(x){
      return(c(rep(1,x),rep(0,n1-x)))
    }))
    t=rep(1:50,each=N)
    m=glm(y1~t,family="binomial")
    b0=coefficients(m)[1]
    b2=coefficients(m)[2]
    return(list(b0=b0,b2=b2,t1=ests_bc$t1,t1t2=ests_bc$t1*invlogit(b0+b2*tps)))
    },mc.cores=cores)
 
  #for(i in 1:its){
  #  ids_star = sample(ids,N,replace=T)
  #  df_star = df
  #  for(j in 1:N){
  #    df_star[(1+(j-1)*t):(j*t),] = df[df$id==ids_star[j],]
  #  }
  #  df_star$id = df$id
  #  ests_star = get_simp_est(df_star)
  #  t1t2[i,]=get_bias_corr_est(ests_star$b0,ests_star$t1,ests_star$b2,t3,tps)
  #}
  t1 = unlist(t1t2_list)[c(TRUE,rep(FALSE,52))]
  b0 = unlist(t1t2_list)[c(FALSE,TRUE,rep(FALSE,51))]
  b2 = unlist(t1t2_list)[c(FALSE,FALSE,TRUE)]
  t1t2 = matrix(unlist(t1t2_list)[c(FALSE,rep(FALSE,50),rep(TRUE,50))],nrow=its,byrow=T)
  return(list(t1=t1,t2=t2,t1t2=t1t2))
}



