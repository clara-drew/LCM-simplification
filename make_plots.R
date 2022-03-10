library(ggplot2)
library(gridExtra)
source("simulation_funs.R")

scenario <- "2b" # Set to given scenario

source(paste("parameters/s",scenario,"-beta0-parameters.R",sep=""))

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

MSE=rep(NA,length(maxT))
CP=rep(NA,length(maxT))
Bias=rep(NA,length(maxT))

MSE_s=rep(NA,length(maxT))
CP_s=rep(NA,length(maxT))
Bias_s=rep(NA,length(maxT))

MSE_BC=rep(NA,length(maxT))
CP_BC=rep(NA,length(maxT))
Bias_BC=rep(NA,length(maxT))

MSE_s_BC=rep(NA,length(maxT))
CP_s_BC=rep(NA,length(maxT))
Bias_s_BC=rep(NA,length(maxT))

MSE_mle=rep(NA,length(maxT))
CP_mle=rep(NA,length(maxT))
Bias_mle=rep(NA,length(maxT))

beta <- c(0,0.005,0.01,0.05)
plotNames <- c("-theta","-alphaFirst","-alphaMiddle","-alphaLast")
betaS = c("0","005","01","05")
mid=rep(2:4,each=2)

grobList=list()

for(k in 1:4){
  
  max_mse <- 0
  min_bias <- 0
  max_bias <- 0

  for(i in 1:4){
  
    for(j in 1:6){
      
      alpha <- alpha_list[[j]]
      
      if(k==1){
        
        df_theta = read.csv(paste("s",scenario,"-beta",betaS[i],"/df_","theta_",maxT[j],".csv",sep=""))
        theta_hat[[j]] = df_theta$theta_hat
        theta_hat_s[[j]] = df_theta$theta_hat_s
        theta_hat_BC[[j]] = df_theta$theta_hat_BC
        theta_hat_s_BC[[j]] = df_theta$theta_hat_s_BC
        theta_hat_MLE[[j]] = df_theta$theta_hat_MLE
        df_theta = read.csv(paste("s",scenario,"-beta",betaS[i],"/df_","theta_cp_",maxT[j],".csv",sep=""))
        theta_hat_CP[[j]] = df_theta$theta_hat
        theta_hat_s_CP[[j]] = df_theta$theta_hat_s
        theta_hat_BC_CP[[j]] = df_theta$theta_hat_BC
        theta_hat_s_BC_CP[[j]] = df_theta$theta_hat_s_BC
        theta_hat_MLE_CP[[j]] = df_theta$theta_hat_MLE
        
        eTheta = theta*mean(1-prod(1-alpha))+(1-theta)*mean(1-(1-beta[i])^maxT[j])
        varTheta = eTheta*(1-eTheta)/100
        
        Bias[j]=eTheta-theta
        CP[j]=mean(theta_hat_CP[[j]])
        MSE[j]=varTheta+Bias[j]^2
        
        MSE_BC[j]=mean((theta_hat_BC[[j]]-theta)^2)
        CP_BC[j]=mean(theta_hat_BC_CP[[j]])
        Bias_BC[j]=mean(theta_hat_BC[[j]])-theta
        
        eTheta=theta*mean(1-prod(1-alpha)-sum(unlist(lapply(2:maxT[j],function(x) alpha[x]*prod(1-alpha[-x]))))) + (1-theta)*mean(1-sum(dbinom(0:1,maxT[j]-1,beta[i])*(1-beta[i])))
        varTheta=eTheta*(1-eTheta)/100
        
        Bias_s[j]=eTheta-theta
        CP_s[j]=mean(theta_hat_s_CP[[j]])
        MSE_s[j]=varTheta+Bias_s[j]^2
        
        MSE_s_BC[j]=mean((theta_hat_s_BC[[j]]-theta)^2)
        CP_s_BC[j]=mean(theta_hat_s_BC_CP[[j]])
        Bias_s_BC[j]=mean(theta_hat_s_BC[[j]])-theta
        
        MSE_mle[j]=mean((theta_hat_MLE[[j]]-theta)^2,na.rm=T)
        CP_mle[j]=mean(theta_hat_MLE_CP[[j]],na.rm=T)
        Bias_mle[j]=mean(theta_hat_MLE[[j]],na.rm=T)-theta
        
      }
      
      else{
        
        ind <- c(1,mid[j],maxT[j])
        
        df_alpha = read.csv(paste("s",scenario,"-beta",betaS[i],"/df_","alpha_",maxT[j],".csv",sep=""))
        alpha_hat[[j]] = subset(df_alpha, select=c(paste("alpha_hat_",1:maxT[j],sep="")))
        alpha_hat_s[[j]] = subset(df_alpha, select=c(paste("alpha_hat_s_",1:maxT[j],sep="")))
        alpha_hat_BC[[j]] = subset(df_alpha, select=c(paste("alpha_hat_BC_",1:maxT[j],sep="")))
        alpha_hat_s_BC[[j]] = subset(df_alpha, select=c(paste("alpha_hat_s_BC_",1:maxT[j],sep="")))
        alpha_hat_MLE[[j]] = subset(df_alpha, select=c(paste("alpha_hat_MLE_",1:maxT[j],sep="")))
        
        df_alpha = read.csv(paste("s",scenario,"-beta",betaS[i],"/df_","alpha_cp_",maxT[j],".csv",sep=""))
        alpha_hat_CP[[j]] = subset(df_alpha, select=c(paste("alpha_hat_",1:maxT[j],sep="")))
        alpha_hat_s_CP[[j]] = subset(df_alpha, select=c(paste("alpha_hat_s_",1:maxT[j],sep="")))
        alpha_hat_BC_CP[[j]] = subset(df_alpha, select=c(paste("alpha_hat_BC_",1:maxT[j],sep="")))
        alpha_hat_s_BC_CP[[j]] = subset(df_alpha, select=c(paste("alpha_hat_s_BC_",1:maxT[j],sep="")))
        alpha_hat_MLE_CP[[j]] = subset(df_alpha, select=c(paste("alpha_hat_MLE_",1:maxT[j],sep="")))
        
        t <- ind[k-1]
        
        eAlpha = (theta*alpha+(1-theta)*beta[i])/(theta*mean(1-prod(1-alpha))+(1-theta)*mean(1-(1-beta[i])^maxT[j]))
        varAlpha = eAlpha[t]*(1-eAlpha[t])/100
        
        Bias[j]=eAlpha[t]-alpha[t]
        CP[j]=mean(alpha_hat_CP[[j]][,t])
        MSE[j]=varAlpha+Bias[j]^2
        
        MSE_BC[j]=mean((alpha_hat_BC[[j]][,t]-alpha[t])^2)
        CP_BC[j]=mean(alpha_hat_BC_CP[[j]][,t])
        Bias_BC[j]=mean(alpha_hat_BC[[j]][,t])-alpha[t]
        
        eTheta=theta*mean(1-prod(1-alpha)-sum(unlist(lapply(2:maxT[j],function(x) alpha[x]*prod(1-alpha[-x]))))) + (1-theta)*mean(1-sum(dbinom(0:1,maxT[j]-1,beta[i])*(1-beta[i])))
        eAlpha = (theta*alpha[t]+(1-theta)*beta[i])/eTheta
        varAlpha = eAlpha*(1-eAlpha)/100
        
        Bias_s[j]=eAlpha-alpha[t]
        CP_s[j]=mean(alpha_hat_s_CP[[j]][,1])
        MSE_s[j]=varAlpha+Bias_s[j]^2
        
        MSE_s_BC[j]=mean((alpha_hat_s_BC[[j]][,t]-alpha[t])^2)
        CP_s_BC[j]=mean(alpha_hat_s_BC_CP[[j]][,t])
        Bias_s_BC[j]=mean(alpha_hat_s_BC[[j]][,t])-alpha[t]
        
        MSE_mle[j]=mean((alpha_hat_MLE[[j]][,t]-alpha[t])^2,na.rm=T)
        CP_mle[j]=mean(alpha_hat_MLE_CP[[j]][,t],na.rm=T)
        Bias_mle[j]=mean(alpha_hat_MLE[[j]][,t],na.rm=T)-alpha[t]
        
      }  
    
    
    }          
  
    df_mse= data.frame(MSE=c(MSE,MSE_BC,MSE_s,MSE_s_BC,MSE_mle),time_BCoints=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
    df_cp= data.frame(CP=c(CP,CP_BC,CP_s,CP_s_BC,CP_mle),time_BCoints=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
    df_bias= data.frame(Bias=c(Bias,Bias_BC,Bias_s,Bias_s_BC,Bias_mle),time_BCoints=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
    max_mse <- max(max_mse, max(df_mse$MSE)+0.01)
    max_bias <- max(max_bias, max(df_bias$Bias)+0.01)
    min_bias <- min(min_bias, min(df_bias$Bias)-0.01)
    
    if(i==4){
    
      p1=ggplot(df_mse,aes(x=time_BCoints,y=MSE,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests") + theme(legend.text=element_text(size=10))+theme(aspect.ratio=4/3)+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))
    
      p2=ggplot(df_cp,aes(x=time_BCoints,y=CP,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests")+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))+ylim(0,1)
    
      p3=ggplot(df_bias,aes(x=time_BCoints,y=Bias,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests")+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))
    
    }
  
    else{
    
      p1=ggplot(df_mse,aes(x=time_BCoints,y=MSE,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
      p2=ggplot(df_cp,aes(x=time_BCoints,y=CP,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+ylim(0,1)+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
      p3=ggplot(df_bias,aes(x=time_BCoints,y=Bias,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
    }
  
    grobList[[1+(i-1)]]=p1
    grobList[[5+(i-1)]]=p2
    grobList[[9+(i-1)]]=p3
  
  
  }
  
  for(i in 1:4){
    grobList[[1+(i-1)]] <- grobList[[1+(i-1)]] + ylim(0,max_mse)
    grobList[[9+(i-1)]] <- grobList[[9+(i-1)]] + ylim(min_bias,max_bias)
  }

  lay = rbind(c(1,1,2,2,3,3,4,4,4),
            c(5,5,6,6,7,7,8,8,8),
            c(9,9,10,10,11,11,12,12,12))

  g=arrangeGrob(grobs=grobList,ncol=4,widths=c(1,1,1,1.5),heights=c(1,1,1,1))

  ggsave(paste("s",scenario,plotNames[k],".png",sep=""),g,width=10,height=10)

}    