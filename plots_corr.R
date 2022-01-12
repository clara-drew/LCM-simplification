library(ggplot2)
#system("unzip t30.zip")
#system("unzip t3005.zip")
#system("unzip t301.zip")
#system("unzip t305.zip")

source("../EM_functions_MAR.R")

maxT = 3:8

t3=c(0,0.005,0.01,0.05)
t1=0.3
b0=2
b2=-(2+1.5)/(maxT-1)

# theta3=0

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

MSE=rep(NA,length(maxT))
CP=rep(NA,length(maxT))
Bias=rep(NA,length(maxT))

MSE_s=rep(NA,length(maxT))
CP_s=rep(NA,length(maxT))
Bias_s=rep(NA,length(maxT))

MSE_p=rep(NA,length(maxT))
CP_p=rep(NA,length(maxT))
Bias_p=rep(NA,length(maxT))

MSE_s_p=rep(NA,length(maxT))
CP_s_p=rep(NA,length(maxT))
Bias_s_p=rep(NA,length(maxT))

MSE_h=rep(NA,length(maxT))
CP_h=rep(NA,length(maxT))
Bias_h=rep(NA,length(maxT))

theta3 = c("0","005","01","05")

grobList=list()

for(i in 1:4){
  
  for(j in 1:6){
    
    t2=invlogit(b0+b2[j]*(0:maxT[j]))
    
    df_t1 = read.csv(paste("ct3",theta3[i],"/df_","t1_",maxT[j],".csv",sep=""))
    df_t2 = read.csv(paste("ct3",theta3[i],"/df_","t2_",maxT[j],".csv",sep=""))
    
    df_t1_fm = read.csv(paste("ct3",theta3[i],"/df_","t1_",maxT[j],"_fm.csv",sep=""))
    df_t2_fm = read.csv(paste("ct3",theta3[i],"/df_","t2_",maxT[j],"_fm.csv",sep=""))
    
    t1_til[[j]] = df_t1$t1_til
    t1_til_s[[j]] = df_t1$t1_til_s
    t2_til[[j]] = subset(df_t2, select=c(paste("t2_til_",1:maxT[j],sep="")))
    t2_til_s[[j]] = subset(df_t2, select=c(paste("t2_til_s_",1:maxT[j],sep="")))
    t1_til_p[[j]] = df_t1$t1_til_p
    t1_til_s_p[[j]] = df_t1$t1_til_s_p
    t2_til_p[[j]] = subset(df_t2, select=c(paste("t2_til_p_",1:maxT[j],sep="")))
    t2_til_s_p[[j]] = subset(df_t2, select=c(paste("t2_til_s_p_",1:maxT[j],sep="")))
    t1_hat[[j]] = df_t1_fm$t1_hat
    t2_hat[[j]] = subset(df_t2_fm, select=c(paste("t2_hat_",1:maxT[j],sep="")))
    
    df_t1 = read.csv(paste("ct3",theta3[i],"/df_","t1_cp_",maxT[j],".csv",sep=""))
    df_t2 = read.csv(paste("ct3",theta3[i],"/df_","t2_cp_",maxT[j],".csv",sep=""))
    
    df_t1_fm = read.csv(paste("ct3",theta3[i],"/df_","t1_cp_",maxT[j],"_fm.csv",sep=""))
    df_t2_fm = read.csv(paste("ct3",theta3[i],"/df_","t2_cp_",maxT[j],"_fm.csv",sep=""))
    
    t1_til_CP[[j]] = df_t1$t1_til
    t1_til_s_CP[[j]] = df_t1$t1_til_s
    t2_til_CP[[j]] = subset(df_t2, select=c(paste("t2_til_",1:maxT[j],sep="")))
    t2_til_s_CP[[j]] = subset(df_t2, select=c(paste("t2_til_s_",1:maxT[j],sep="")))
    t1_til_p_CP[[j]] = df_t1$t1_til_p
    t1_til_s_p_CP[[j]] = df_t1$t1_til_s_p
    t2_til_p_CP[[j]] = subset(df_t2, select=c(paste("t2_til_p_",1:maxT[j],sep="")))
    t2_til_s_p_CP[[j]] = subset(df_t2, select=c(paste("t2_til_s_p_",1:maxT[j],sep="")))
    t1_hat_CP[[j]] = df_t1_fm$t1_hat
    t2_hat_CP[[j]] = subset(df_t2_fm, select=c(paste("t2_hat_",1:maxT[j],sep="")))
    
    Et1 = t1*mean(1-prod(1-t2))+(1-t1)*mean(1-(1-t3[i])^maxT[j])
    Vart1 = Et1*(1-Et1)/100
    
    #MSE[j]=mean((t1_til[[j]]-t1)^2)
    CP[j]=mean(t1_til_CP[[j]])
    #Bias[j]=mean(t1_til[[j]])-t1
    
    Bias[j]=Et1-t1
    MSE[j]=Vart1+Bias[j]^2
    
    MSE_p[j]=mean((t1_til_p[[j]]-t1)^2)
    CP_p[j]=mean(t1_til_p_CP[[j]])
    Bias_p[j]=mean(t1_til_p[[j]])-t1
    
    Et1=t1*mean(1-prod(1-t2)-sum(unlist(lapply(2:maxT[j],function(x) t2[x]*prod(1-t2[-x]))))) + (1-t1)*mean(1-sum(dbinom(0:1,maxT[j]-1,t3[i])*(1-t3[i])))
    Vart1=Et1*(1-Et1)/100
    
    #MSE_s[j]=mean((t1_til_s[[j]]-t1)^2)
    CP_s[j]=mean(t1_til_s_CP[[j]])
    #Bias_s[j]=mean(t1_til_s[[j]])-t1
    
    Bias_s[j]=Et1-t1
    MSE_s[j]=Vart1+Bias_s[j]^2
    
    MSE_s_p[j]=mean((t1_til_s_p[[j]]-t1)^2)
    CP_s_p[j]=mean(t1_til_s_p_CP[[j]])
    Bias_s_p[j]=mean(t1_til_s_p[[j]])-t1
    
    MSE_h[j]=mean((t1_hat[[j]]-t1)^2,na.rm=T)
    CP_h[j]=mean(t1_hat_CP[[j]],na.rm=T)
    Bias_h[j]=mean(t1_hat[[j]],na.rm=T)-t1
    
  }          
  
  df_mse= data.frame(MSE=c(MSE,MSE_p,MSE_s,MSE_s_p,MSE_h),time_points=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
  df_cp= data.frame(CP=c(CP,CP_p,CP_s,CP_s_p,CP_h),time_points=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
  df_bias= data.frame(Bias=c(Bias,Bias_p,Bias_s,Bias_s_p,Bias_h),time_points=rep(maxT,5),Estimate=c(rep("z",2*length(maxT)),rep("z*",2*length(maxT)),rep("MLE",length(maxT))),bc=as.factor(c(rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)),rep(1,length(maxT)),rep(0,length(maxT)))))
  
  if(i==4){
    
    p1=ggplot(df_mse,aes(x=time_points,y=MSE,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests")+ylim(0,.1) + theme(legend.text=element_text(size=10))+theme(aspect.ratio=4/3)+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))
    
    p2=ggplot(df_cp,aes(x=time_points,y=CP,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests")+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))+ylim(0,1)
    
    p3=ggplot(df_bias,aes(x=time_points,y=Bias,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+scale_linetype_discrete(name="Bias Corrected",breaks=c(0,1),labels=c("No","Yes"))+xlab("Number of tests")+ylim(-.1,.3)+theme(legend.text=element_text(size=8),legend.title=element_text(size=10))
    
    #ggsave(paste("mse_t1_",theta3[i],"_fix.png",sep=""), p1, height=2.78,width=4,units="in")
    #ggsave(paste("cp_t1_",theta3[i],"_fix.png",sep=""), p2, height=2.78,width=4,units="in")
    #ggsave(paste("bias_t1_",theta3[i],"_fix.png",sep=""), p3, height=2.78,width=4,units="in")
    
  }
  
  else{
    
    p1=ggplot(df_mse,aes(x=time_points,y=MSE,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+ylim(0,.1)+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
    p2=ggplot(df_cp,aes(x=time_points,y=CP,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+ylim(0,1)+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
    p3=ggplot(df_bias,aes(x=time_points,y=Bias,color=Estimate,linetype=as.factor(bc))) + geom_line(alpha=0.8)+xlab("Number of tests")+ylim(-.1,.3)+theme(legend.position="none")+theme(aspect.ratio=4/3)
    
    #ggsave(paste("mse_t1_",theta3[i],"_fix.png",sep=""), p1, height=2.78,width=2.7,units="in")
    #ggsave(paste("cp_t1_",theta3[i],"_fix.png",sep=""), p2, height=2.78,width=2.7,units="in")
    #ggsave(paste("bias_t1_",theta3[i],"_fix.png",sep=""), p3, height=2.78,width=2.7,units="in")
    
    
  }
  
  grobList[[1+(i-1)]]=p1
  grobList[[5+(i-1)]]=p2
  grobList[[9+(i-1)]]=p3
  
  
}

lay = rbind(c(1,1,2,2,3,3,4,4,4),
            c(5,5,6,6,7,7,8,8,8),
            c(9,9,10,10,11,11,12,12,12))

g=arrangeGrob(grobs=grobList,ncol=4,widths=c(1,1,1,1.5),heights=c(1,1,1,1))

dev.off()


ggsave("plots_t1_corr_fm.png",g,width=10,height=10)

middle=rep(2:4,each=2)