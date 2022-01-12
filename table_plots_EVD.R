semen = read.csv("semen_cep1.csv")
source("EM_functions_MAR.R")
library(dplyr)
library(xtable)
library(ggplot2)
library(parallel)
library(mgcv)
library(lessR)
library(flexmix)
library(npreg)

lambda <- 0.3

invlogit <- function(x) exp(x)/(1+exp(x))

semen <- read.csv("semen_cep1.csv")
semen <- semen[complete.cases(semen),]
semen$id <- as.factor(semen$id)

pids <- unique(semen$id)
n <- length(pids)
maxT=max(semen$t)

semen_short <- semen %>% group_by(id) %>% filter(row_number() == 1)
semen$z <- NA

k=1
for(i in 1:dim(semen_short)[1]){
  times<-as.numeric(table(semen$id)[i])
  semen$z[k:(k+times-1)] <- rep(semen_short$y[i]+1,times)
  k<-k+times
}

m1 <- flexmix(cbind(y,1-y) ~ 1 | id, data = semen, k=2, cluster=semen$z,
              model = FLXMRglmfix(nested = list(k=c(1,1), formula = c(~0,~t)),family="binomial"))


refit_m1 <- refit(m1)

coefs <- refit_m1@coef
Cov <- refit_m1@vcov

b0 <- coefs[3]
b1 <- coefs[4]
b2 <- coefs[1]
b3 <- coefs[2]

tps <- c(1,11,21,31)
tps2 <- c(9,17,32,41)

t1 <- invlogit(b1)
t2 <- invlogit(b0+b2*tps)
t3 <- invlogit(b3)

t2_EM = invlogit(b0+b2*(1:50))

EM = c(t1,b0,b2,t2,t3)

colnames(Cov) <- c("b2","b3","b0","b1")
rownames(Cov) <- colnames(Cov)

Cov <- Cov[c(3,4,1,2),c(3,4,1,2)]

y <- semen$y
t <- semen$t
id <- semen$id

lowerExact = function(p,n){
  obs = !is.na(p)
  p = p[obs]
  n = n[obs]
  x = p*n
  len = length(x)
  lower = rep(NA,len)
  for(i in 1:len){
    lower[i]=binom.test(x[i],n[i])$conf.int[1]
  }
  return(lower)
}

upperExact = function(p,n){
  obs = !is.na(p)
  p = p[obs]
  n = n[obs]
  x = p*n
  len = length(x)
  upper = rep(NA,len)
  for(i in 1:len){
    upper[i]=binom.test(x[i],n[i])$conf.int[2]
  }
  return(upper)
}


dqdb = function(t) return(c(invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),0,t*invlogit(b0+b2*t)*(1-invlogit(b0+b2*t)),0))

dgdb = function() return(c(0,0,0,invlogit(b3)*(1-invlogit(b3))))

se_t1 = sqrt(c(0,invlogit(b1)*(1-invlogit(b1)),0,0) %*% Cov %*% c(0,invlogit(b1)*(1-invlogit(b1)),0,0))
se_t2 = unlist(lapply(1:50, function(x) sqrt(dqdb(x) %*% Cov %*% dqdb(x))))
se_t3 = sqrt(c(0,0,0,invlogit(b3)*(1-invlogit(b3))) %*% Cov %*%  c(0,0,0,invlogit(b3)*(1-invlogit(b3))))

SE = c(se_t1,sqrt(Cov[1,1]),sqrt(Cov[3,3]),se_t2[tps],se_t3)

lower = EM - qnorm(.975)*SE

upper = EM + qnorm(.975)*SE

EM_row = paste(round(EM,3)," (",round(lower,3),", ",round(upper,3),")",sep="")

lowert2 = invlogit(b0+b2*(1:50))-qnorm(.975)*se_t2
uppert2 = invlogit(b0+b2*(1:50))+qnorm(.975)*se_t2

lowerMLE <- lowert2[tps2]
upperMLE <- uppert2[tps2]

N_list = list()

# Zi - simple

df <- data.frame(t,id,y)
z = ifelse(aggregate(df$y,by=list(id=df$id),FUN=sum,na.rm=T)[,2]>0,1,0)
df_zi = left_join(data.frame(id=unique(df$id),z=z),df)

n2=left_join(data.frame(t=1:50),aggregate(df_zi[!is.na(df_zi$y),]$z,by=list(t=df_zi[!is.na(df_zi$y),]$t),FUN=sum))[,2]

N_list[[1]]=n2[!is.na(n2) & n2!=0]

N = unlist(c(n,n2[tps]))

t1_samp=mean(z)
t2_samp = left_join(data.frame(t=1:50),aggregate(y~t,data=df_zi[df_zi$z==1,],FUN=mean,na.rm=T))[,2]

Zi = c(t1_samp, t2_samp[tps])

#SE = c(sqrt(Zi[-(length(Zi))]*(1-Zi[-(length(Zi))])/N),NA)

#lower = Zi - qnorm(.975)*SE

#upper = Zi + qnorm(.975)*SE

lower = lowerExact(Zi,N)

upper = upperExact(Zi,N)

#lowert2 = c(lowert2,lowerExact(t2_samp,n2))
#uppert2 = c(uppert2,upperExact(t2_samp,n2))

Zi_row = paste(round(Zi,3)," (",round(lower,3),", ",round(upper,3),")",sep="")

Zi_row = c(Zi_row[1],"","",Zi_row[2:5],0)

# Zi* - simple

df_obs = df[!is.na(df$y),]
df_obs$id <- as.character(df_obs$id)
z1 = df_obs[!duplicated(df_obs$id),]$y
z2 = ifelse(aggregate(df$y,by=list(id=df$id),FUN=sum,na.rm=T)[,2]>1,1,0)
z = ifelse(z1+z2>0,1,0)
df_zi = left_join(data.frame(id=unique(df$id),z=z),df)

n2=left_join(data.frame(t=1:50),aggregate(df_zi[!is.na(df_zi$y),]$z,by=list(t=df_zi[!is.na(df_zi$y),]$t),FUN=sum))[,2]

N_list[[2]]=n2[!is.na(n2) & n2!=0]

N = c(n,n2[tps],dim(df_zi[df_zi$z==0 & !is.na(df_zi$y),])[1])

t1_samp_s= mean(z)
t2_samp_s= left_join(data.frame(t=1:50),aggregate(y~t,data=df_zi[df_zi$z==1,],FUN=mean,na.rm=T))[,2]
t3_samp_s = mean(df_zi[df_zi$z==0,]$y,na.rm=T)

dfy = df_zi[df_zi$z==1 & !is.na(df_zi$y),]

m = gam(dfy$y~s(dfy$t))

dfy$fit = m$fitted.values

se = c(se_t2,predict(m,se.fit=T)$se.fit)

smoothdf = data.frame(fit=c(invlogit(b0+b2*(1:50)),dfy$fit),t=c(1:50,dfy$t),Estimate=c(rep("MLE",50),rep("z*",dim(dfy)[1])))

p=ggplot(smoothdf,aes(t,fit)) + geom_ribbon(aes(ymin=fit-qnorm(.975)*se,ymax=fit+qnorm(.975)*se,fill=Estimate),alpha=0.3) + geom_point(data=data.frame(t=c(1:50),fit=t2_samp_s,Estimate=rep("z*",50)),aes(t,fit,col=Estimate)) + geom_line(aes(col=Estimate),size=1)+xlab("Time")+ylab("alpha")+geom_jitter(data=data.frame(t=df_obs$t,fit=rep(-.2,dim(df_obs)[1])),aes(t,fit),alpha=0.1,width=0.3,height=0.03)



Zi = c(t1_samp_s, t2_samp_s[tps],t3_samp_s)


lower = lowerExact(Zi,N)

upper = upperExact(Zi,N)

lowerZs <- lowerExact(t2_samp_s[tps2],n2[tps2])
upperZs <- upperExact(t2_samp_s[tps2],n2[tps2])


Zis_row = paste(round(Zi,3)," (",round(lower,3),", ",round(upper,3),")",sep="")

Zis_row = c(Zis_row[1],"","",Zis_row[2:6])



# Bias-corrected

gp = rep(NA,n)

for(i in 1:n){
  tp = df_obs[df_obs$id==pids[i],]$t
  gp[i]=prod(1-t2_samp[tp],na.rm=T)
}

Zi1 = c(t1_samp/mean(1-gp),NA,NA,(t2_samp*mean(1-gp))[tps],0)

t2_samp_p = t2_samp*mean(1-gp)

pz10 = rep(NA,n)
pz11 = rep(NA,n)
pz111 = matrix(nrow=n,ncol=maxT)
pz110 = rep(NA,n)

for(i in 1:n){
  tp = df_obs[df_obs$id==pids[i],]$t
  Ti = length(tp)
  pz10[i] = 1-sum(dbinom(0:1,Ti-1,t3_samp_s)*(1-t3_samp_s))
  pz11[i] = 1-prod(1-t2_samp_s[tp],na.rm=T)-sum(unlist(lapply(tp[-1],FUN=function(x) t2_samp_s[x]*prod(1-t2_samp_s[-x],na.rm=T))),na.rm=T)
  pz111[i,tp] = c(1,unlist(lapply(tp[-1],FUN = function(x) 1-prod(1-t2_samp_s[tp[-x]]))))
  pz110[i] = 1-(1-t3_samp_s)^(Ti-1)
}

t1tilp=(t1_samp_s-mean(pz10))/mean(pz11-pz10)
t2tilp=(t2_samp_s*(mean(pz11)*t1tilp+mean(pz10)*(1-t1tilp))-
        mean(pz110)*t3_samp_s*t1tilp)/(colMeans(pz111,na.rm=T)*t1tilp)
t3tilp=(t3_samp_s*(mean(1-pz11)*t1tilp+mean(1-pz10)*(1-t1tilp))-
        mean(colMeans(1-pz111,na.rm=T)*t2tilp*(1-t1tilp),na.rm=T))/mean((1-pz110)*t1tilp)

t2_samp_s_p = t2tilp

Zi2 = c(t1tilp,NA,NA,t2tilp[tps],t3tilp)


iter = function(x){

  ids = sample(pids,n,replace=T)
  dfB=data.frame(t=NULL,id=NULL,y=NULL)
  for(i in 1:n){
    nl = dim(dfB)[1]
    dfB = rbind(dfB,df_obs[df_obs$id==ids[i],])
    ni = dim(df_obs[df_obs$id==ids[i],])[1]
    dfB$id[(nl+1):(nl+ni)]=as.character(i)
  }
  z_samp = ifelse(aggregate(dfB$y,by=list(id=dfB$id),FUN=sum)[,2]>0,1,0)
  t11 = mean(z_samp)
  id = 1:n
  df_boot = left_join(data.frame(id=as.character(id),z1=z_samp),dfB)
  z_samp1 = dfB[!duplicated(dfB$id),]$y
  z_samp2 = ifelse(aggregate(dfB$y,by=list(id=dfB$id),FUN=sum)[,2]>1,1,0)
  z_samp = ifelse(z_samp1+z_samp2>0,1,0)
  t12 = mean(z_samp)
  df_boot = left_join(data.frame(id=as.character(1:n),z2=z_samp),df_boot)

  t21 = left_join(data.frame(t=1:50),aggregate(df_boot[df_boot$z1==1,]$y,by=list(t=df_boot[df_boot$z1==1,]$t),FUN=mean))[,2]
  t22 = left_join(data.frame(t=1:50),aggregate(df_boot[df_boot$z2==1,]$y,by=list(t=df_boot[df_boot$z2==1,]$t),FUN=mean))[,2]

  t3 = mean(df_boot[df_boot$z2==0,]$y)

  gp = rep(NA,length(id))
  pz10 = rep(NA,n)
  pz11 = rep(NA,n)
  pz111 = matrix(nrow=n,ncol=maxT)
  pz110 = rep(NA,n)

  for(i in 1:n){
    tp = df_boot[df_boot$id==i,]$t
    Ti = length(tp)
    pz10[i] = (1-sum(dbinom(0:1,Ti-1,t3)*(1-t3)))
    pz11[i] = (1-prod(1-t22[tp],na.rm=T)-sum(unlist(lapply(tp[-1],FUN=function(x) t22[x]*prod(1-t22[-x],na.rm=T))),na.rm=T))       
    pz111[i,tp] = c(1,unlist(lapply(tp[-1],FUN = function(x) 1-prod(1-t22[tp[-x]]))))
    pz110[i] = 1-(1-t3)^(Ti-1)
    gp[i]=prod(1-t21[tp],na.rm=T)
  }

  t11p = t11/mean(1-gp)
  t21p = t21*mean(1-gp)
  t12p=(t12-mean(pz10))/mean(pz11-pz10)
  t22p=(t22*(mean(pz11)*t12p+mean(pz10)*(1-t12p))-
        mean(pz110)*t3*t12p)/(colMeans(pz111,na.rm=T)*t12p)
  t3p=(t3*(mean(1-pz11)*t12p+mean(1-pz10)*(1-t12p))-
        mean(colMeans(1-pz111,na.rm=T)*t22*(1-t12p),na.rm=T))/mean((1-pz110)*t12p)
  obs=!is.na(t22p)
  
  n2=left_join(data.frame(t=1:50),aggregate(df_zi[!is.na(df_zi$y),]$z,by=list(t=df_zi[!is.na(df_zi$y),]$t),FUN=sum))[,2]
  
  t22s=npreg::ss(y=t22p[obs],x=c(1:50)[obs],w=n2[obs],spar=lambda)$y
  #t22s=gam(t22p[obs]~s((1:50)[obs]))$fitted.values
  t22s=left_join(data.frame(t=1:50),data.frame(t=(1:50)[obs],t22s=t22s))$t22s

  return(list(t11p=t11p,t21p=t21p,t12p=t12p,t22p=t22p,t3p=t3p,t22s=t22s))

}

bootCI = function(its=1000,cores=3){
  est_list = mclapply(1:its,FUN=iter,mc.cores=cores)
  N = length(unlist(est_list))/its-50
  CIs = list()
  SEs = list()
  ind = rep(FALSE,N+50)
  for(i in 1:N){
    indi = ind
    indi[i] = TRUE
    CIs[[i]] = quantile(unlist(est_list)[indi],probs=c(0.025,0.975),na.rm=T)
  }
  for(i in (N+1):(N+50)){
   indi = ind
   indi[i] = TRUE
   SEs[[i]] = sd(unlist(est_list)[indi],na.rm=T)
  }
  return(list(CIs=CIs,SEs=SEs))
}

boot = bootCI(its=100,cores=1)

boot_ci = boot$CIs
boot_se = unlist(boot$SEs)

lower_full= unlist(boot_ci)[c(TRUE,FALSE)]

upper_full= unlist(boot_ci)[c(FALSE,TRUE)]

ind1 = c(1,tps+1)
ind2 = c(52,(c(tps,51)+52))

lower = c(lower_full[1],NA,NA,lower_full[tps+1],lower_full[52],NA,NA,lower_full[c(tps,51)+52])

upper = c(upper_full[1],NA,NA,upper_full[tps+1],upper_full[52],NA,NA,upper_full[c(tps,51)+52])


#lowert2 = c(lowert2,lower_full[c(2:51,53:102)])
#uppert2 = c(uppert2,upper_full[c(2:51,53:102)])

Zi = c(Zi1[-8],Zi2)

both_row = paste(round(Zi,3)," (",round(lower,3),", ",round(upper,3),")",sep="")


ind1 = c(1,tps2+1)
ind2 = c(52,(c(tps2,51)+52))

lower = c(lower_full[1],NA,NA,lower_full[tps+1],lower_full[52],NA,NA,lower_full[c(tps,51)+52])

upper = c(upper_full[1],NA,NA,upper_full[tps+1],upper_full[52],NA,NA,upper_full[c(tps,51)+52])

lowerZsp <- lower[11:14]
upperZsp <- upper[11:14]


Zip_row = c(both_row[1],"","",both_row[4:7],"")

Zisp_row = c(both_row[8],"","",both_row[11:15])

table = rbind(EM_row,Zi_row,Zis_row,Zip_row,Zisp_row)

table = table[,-c(2:3)]

colnames(table) = c("$\\theta$","$\\alpha_1$","$\\alpha_{11}$","$\\alpha_{21}$","$\\alpha_{31}$","$\\beta$")

rownames(table) = c("MLE","$z_i$","$z_i^*$","$z_i$-BC","$z_i^*$-BC")

print(xtable(table,label="treview",digits=3,align="lcccccc"),table.placement="H",file="prevail_vals.tex",booktabs=T,sanitize.text.function=function(x){x})



t=c(1,5:50)

t2_samp = t2_samp[!is.na(t2_samp)]
t2_samp_s= t2_samp_s[!is.na(t2_samp_s)]
t2_samp_p= t2_samp_p[!is.na(t2_samp_p)]
t2_samp_s_p= t2_samp_s_p[!is.na(t2_samp_s_p)]

m1 = glm(t2_samp~t,family="binomial",weights=N_list[[1]])
m2 = glm(t2_samp_s~t,family="binomial",weights=N_list[[2]])
m3 = glm(t2_samp_p~t,family="binomial",weights=N_list[[1]])

t2_samp_s_p=ifelse(t2_samp_s_p<0,0,ifelse(t2_samp_s_p>1,1,t2_samp_s_p))

m4 = glm(t2_samp_s_p~t,family="binomial",weights=N_list[[2]])

nd = data.frame(t=1:50)

preds = list(predict(m1,newdata=nd,type="response",se.fit=T),predict(m2,newdata=nd,type="response",se.fit=T),predict(m3,newdata=nd,type="response",se.fit=T),predict(m4,newdata=nd,type="response",se.fit=T))

df_points = data.frame(prob=c(t2_samp_s,t2_samp_s_p),t=rep(c(1,5:50),2),c(rep("No",47),rep("Yes",47)),Estimate=rep("z*",94))

df = data.frame(prob=c(preds[[1]]$fit,preds[[2]]$fit,preds[[3]]$fit,preds[[4]]$fit,t2_EM),lower=c(unlist(lapply(preds, FUN = function(x) x$fit-qnorm(.975)*x$se.fit)),lowert2),upper=c(unlist(lapply(preds, FUN = function(x) x$fit+qnorm(.975)*x$se.fit)),uppert2),Estimate=c(rep("z",50),rep("z*",50),rep("z",50),rep("z*",50),rep("MLE",50)),BC=c(rep("No",100),rep("Yes",100),rep("No",50)),t=rep(1:50,5))

p=ggplot(df[df$Estimate!="z",],aes(x=t,y=prob,color=Estimate,linetype=BC)) + geom_line() + xlab("Time") + ylab("Probability") + geom_ribbon(aes(ymin=lower,ymax=upper,x=t,fill=Estimate),alpha=.1,color=NA) + labs(linetype="Bias-Corrected")

#ggsave("glm_alpha_se.png",p)

p=ggplot(df,aes(x=t,y=prob,color=Estimate,linetype=BC)) + geom_line() + xlab("Time") + ylab("Probability") + labs(linetype="Bias-Corrected")

#ggsave("glm_alpha.png",p)


# Spline plots

smoothT2 <- rep(NA,50)

smoothT2[c(1,5:50)] <- npreg::ss(y=t2_samp_s_p,x=c(1,5:50),w=n2[-c(2:4)],spar=lambda)$y

lowerZsp <- (smoothT2-qnorm(0.975)*boot_se)[tps2]
upperZsp <- (smoothT2+qnorm(0.975)*boot_se)[tps2]

smoothT2 <- npreg::ss(y=df_zi[df_zi$z==1,]$y,x=df_zi[df_zi$z==1,]$t,spar=lambda)

Ses <- predict(smoothT2,x=tps2,se.fit=TRUE)[,3]

lowerZs <- predict(smoothT2,x=tps2)[,2]-qnorm(0.975)*Ses
upperZs <- predict(smoothT2,x=tps2)[,2]+qnorm(0.975)*Ses

lower <- ifelse(c(lowerZs, lowerZsp, lowerMLE)<0,0,c(lowerZs, lowerZsp, lowerMLE))
upper <- c(upperZs, upperZsp, upperMLE)

pd <- position_dodge(width=3.5)

df <- data.frame(prob=c(npreg::ss(y=df_zi[df_zi$z==1,]$y,x=df_zi[df_zi$z==1,]$t,spar=lambda)$y, 
          npreg::ss(y=t2_samp_s_p,x=c(1,5:50),w=n2[-c(2:4)],spar=lambda)$y, t2_EM),
          Estimate=c(rep("z*",47*2),rep("MLE",50)),BC=c(rep("No",47),rep("Yes",47),rep("No",50)),
          Time=c(1,5:50,1,5:50,1:50)) 

df_error <- data.frame(Time=rep(tps2,3),lower=lower,upper=upper,prob=df$prob[df$Time %in% tps2],
                       Estimate=df$Estimate[df$Time %in% tps2], BC=df$BC[df$Time %in% tps2])

p=ggplot(df,aes(x=Time,y=prob,color=Estimate,linetype=BC)) + geom_line() + xlab("Time") + 
  ylab("Probability") + labs(linetype="Bias-Corrected") +
 geom_errorbar(data=df_error,aes(ymin=lower,ymax=upper), position=pd)

ggsave("cubic-spline-error.png",p)
