library(ggplot2)
source("../../EM_functions_MAR.R")

maxT = 8
t3 = 0.01
t1=0.3
b0=2
b2 = -(2+1.5)/(maxT-1)
n=c(100,200,300,400,500)
t2 = invlogit(b0+b2*(0:maxT))

df_cp = data.frame(theta=rep(NA,25),SS=rep(n,each=5),Estimate=rep(c("z","z*","z-BC","z*-BC","MLE"),5))

df_me = data.frame(theta=rep(NA,25),SS=rep(n,each=5),Estimate=rep(c("z","z*","z-BC","z*-BC","MLE"),5))

ind = split(1:25,rep(1:5,each=5))

for(i in 1:length(n)){
  
  theta_cp=read.csv(paste("df_t1_cp_",n[i],".csv",sep=""))
  theta_me=read.csv(paste("df_t1_ME_",n[i],".csv",sep=""))
  
  theta_cp_fm=read.csv(paste("df_t1_cp_fm",n[i],".csv",sep=""))
  theta_me_fm=read.csv(paste("df_t1_ME_fm",n[i],".csv",sep=""))

  df_cp$theta[ind[[i]][1]] = mean(theta_cp$t1_til)
  df_cp$theta[ind[[i]][2]] = mean(theta_cp$t1_til_s)
  df_cp$theta[ind[[i]][3]] = mean(theta_cp$t1_til_p)
  df_cp$theta[ind[[i]][4]] = mean(theta_cp$t1_til_s_p)
  df_cp$theta[ind[[i]][5]] = mean(theta_cp_fm$t1_hat)

  df_me$theta[ind[[i]][1]] = mean(theta_me$t1_til)
  df_me$theta[ind[[i]][2]] = mean(theta_me$t1_til_s)
  df_me$theta[ind[[i]][3]] = mean(theta_me$t1_til_p)
  df_me$theta[ind[[i]][4]] = mean(theta_me$t1_til_s_p)
  df_me$theta[ind[[i]][5]] = mean(theta_me_fm$t1_hat)
}

MEs = sqrt((qnorm(.975)^2*.5^2)/n)

df_ref = data.frame(theta=MEs,SS=n)

p1 = ggplot(df_me,aes(x=SS,y=theta,color=Estimate)) + geom_line() + ylab("Margin of Error") + xlab("Sample Size") + geom_line(data=df_ref,aes(x=SS,y=theta),linetype="dashed",color="Black")+theme(legend.position="none")+theme(aspect.ratio=4/3)
p2 = ggplot(df_cp,aes(x=SS,y=theta,color=Estimate)) + geom_line() + ylab("Coverage Probability") + xlab("Sample Size") + geom_hline(yintercept=0.95,color="Black",linetype="dashed")+theme(aspect.ratio=4/3)

g = arrangeGrob(grobs=list(p1,p2),ncol=2,widths=c(1,1.5))

ggsave("theta_SS_test.png",g,width=10,height=5)
