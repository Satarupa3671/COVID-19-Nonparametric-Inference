
### ggplot of CI from sampling_agg_sim.RData


rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)

### set 1
######
Time <- 100
true.para <- list(gamma.true = 0.03,
                  rhoA.true = 0.06,#0.02 (relative smaller NIR but larger than 0.05)
                  rhoH.t.true = c(rep(0.03,20),rep(0.04,40),rep(0.05,40)), 
                  alpha.true = 0.2,
                  delta.t.true = c(rep(0.007,50),rep(0.005,50)),
                  phi.t.true = 0.001+1:Time*0.00005)
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/set1_wd_30.RData")
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/sampling_CI/set1_wd_30/set1_wd_30_sampling_agg_sim.RData")
#########

### set 2
######
Time <- 100

true.para <- list(gamma.true = 0.01,
                  rhoA.true = 0.04,#0.02 (relative smaller NIR but larger than 0.05)
                  rhoH.t.true = c(rep(0.02,20),rep(0.03,40),rep(0.04,40)), 
                  alpha.true = 0.2,
                  delta.t.true = c(rep(0.005,50),rep(0.008,50)),
                  phi.t.true = 0.002+1:Time*0.00005)
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/set2_wd_30.RData")
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/sampling_CI/set2_wd_30/set2_wd_30_sampling_agg_sim.RData")
###########

## set 3:
######

Time <- 100

true.para <- list(gamma.true = 0.01,
                  rhoA.true = 0.02,
                  rhoH.t.true = c(rep(0.02,20),rep(0.05,50),rep(0.04,30)), 
                  alpha.true = 0.18,
                  delta.t.true = c(rep(0.005,30),rep(0.009,40),rep(0.007,30)),
                  phi.t.true = 0.0001+1:Time*0.00005)
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/set3_wd_30.RData")
load("/Users/shutingliao/Documents/UCdavis/Project/epidemic_mech/Dropbox/simulation/aggregation/sampling_CI/set3_wd_30/set3_wd_30_sampling_agg_sim.RData")

###########

sim.BS.N = 1050
level  = .05

s <- list()
## collect gamma

pdf("hist_gamma.pdf")
df_gamma = as.data.frame(cbind(gamma = est.bs$gamma.bs, gamma.true = true.para$gamma.true,
                               gamma.mean = mean(est.bs$gamma.bs)))
s[[1]] <- ggplot(df_gamma, aes(x=gamma,y = ..density..))+
  geom_histogram(color="darkorange", fill="orange")+
  geom_density(color="darkred",adjust=3)+ #set2:adjust=3
  labs(#title = expression(paste('Sampling distribution of ',gamma)), 
       x = expression(gamma), y = "Density") +theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_vline(aes(xintercept = gamma.true),col='blue',size=1)+
  geom_vline(aes(xintercept = gamma.mean),col='red',size=1)
dev.off()

sd(est.bs$gamma.bs)
mean(est.bs$gamma.bs)
sort(est.bs$gamma.bs)[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]





## collect rhoA
df_rhoA = as.data.frame(cbind(rhoA = est.bs$rhoA.bs, rhoA.true = true.para$rhoA.true,
                              rhoA.mean  = mean(est.bs$rhoA.bs)))
pdf('hist_rhoA.pdf')
s[[2]] <- ggplot(df_rhoA, aes(x=rhoA,y = ..density..))+
  geom_histogram(color="darkorange", fill="orange")+
  geom_density(color="darkred")+# set1:,adjust=2
  labs(#title = expression(paste('Sampling distribution of ',rho[A])), 
       x = expression(rho[A]), y = "Density") +theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_vline(aes(xintercept = rhoA.true),col='blue',size=1)+
  geom_vline(aes(xintercept = rhoA.mean),col='red',size=1)
dev.off()

sd(est.bs$rhoA.bs)
mean(est.bs$rhoA.bs)
sort(est.bs$rhoA.bs)[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]


## collect At
At.bs = est.bs$At.bs
At.ci.bs <- t(sapply(1:dim(At.bs)[1], function(i){
  sort(At.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))
At.true <- c(sapply(1:10,function(i){dat.sim[[i]]$data$At[1:(dim(At.bs)[1]/10)]}))

df_At = as.data.frame(cbind(Time = rep(1:(dim(At.bs)[1]/10),10), m.lower = At.ci.bs[,1], 
                            m.upper = At.ci.bs[,2], mu = rowMeans(est.bs$At.bs),
                            true = At.true))

pdf("ci_At.pdf")
p <- list()
for(i in 1:10){
  
  df_At_sub <- df_At[((i-1)*(dim(At.bs)[1]/10)+1):(dim(At.bs)[1]/10*i),]
  p[[i]] <- ggplot(df_At_sub, aes(Time)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu), linetype = 'solid', size = 1,col='red') +
    geom_line(aes(y = true), linetype = 'solid', size = 1,col='blue') +
    labs(#title = bquote('CI of A(t) for state'~.(i)),
      x = "Time", y = as.expression(bquote(paste("A"[t],"-process"~.(i))))) +
    theme_bw()+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
             p[[10]],nrow=4,ncol=3)
dev.off()

## collect Ct
Ct.bs <- est.bs$Ct.bs
Ct.ci.bs <- t(sapply(1:dim(Ct.bs)[1], function(i){
  sort(Ct.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))
Ct.true <- c(sapply(1:10,function(i){dat.sim[[i]]$data$Ct[1:(dim(Ct.bs)[1]/10)]}))

df_Ct = as.data.frame(cbind(Time = rep(1:(dim(Ct.bs)[1]/10),10), m.lower = Ct.ci.bs[,1], 
                            m.upper = Ct.ci.bs[,2], mu = rowMeans(est.bs$Ct.bs),
                            true = Ct.true))

pdf("ci_Ct.pdf")
p <- list()
for(i in 1:10){
  
  df_Ct_sub <- df_Ct[((i-1)*(dim(Ct.bs)[1]/10)+1):(dim(Ct.bs)[1]/10*i),]
  p[[i]] <- ggplot(df_Ct_sub, aes(Time)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu), linetype = 'solid', size = 1) +
    geom_line(aes(y = true), linetype = 'solid', size = 1,col='red') +
    labs(#title = bquote('CI of C(t) for state'~.(i)),
         x = "Time", y = as.expression(bquote(paste("C"[t],"-process"~.(i)))))+
    theme_bw()+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
             p[[10]],nrow=4,ncol=3)
dev.off()


## point-wise CI for rhoH.t, phi.t
rhoH.smooth.bs <- est.bs$rhoH.smooth.bs
#range(rhoH.bs)
rhoH.ci.bs <- t(sapply(1:dim(rhoH.smooth.bs)[1], function(i){
  sort(rhoH.smooth.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))



df_rhoH = as.data.frame(cbind(Time = 1:dim(rhoH.ci.bs)[1], 
                              m.lower = rhoH.ci.bs[,1], 
                              m.upper = rhoH.ci.bs[,2],
                              mu = rowMeans(rhoH.ci.bs),
                              true = true.para$rhoH.t.true[1:dim(rhoH.ci.bs)[1]]))

pdf("ci_rhoH.pdf")
s[[3]] <- ggplot(df_rhoH, aes(Time)) + 
  geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
              alpha=0.7)+
  geom_line(aes(y = mu), linetype = 'solid', size = 1,col='red') +
  geom_line(aes(y = true), linetype = 'solid', size = 1,col='blue') +
  labs(#title =expression(paste("'Pointwise Confidence Interval for ",rho[H],'(t)')),
       x = "Time", y = expression(paste(rho[H],'(t)')))  +theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))
dev.off()


##
phi.smooth.bs <- est.bs$phi.smooth.bs
#range(phi.smooth.bs)
phi.ci.bs <- t(sapply(1:dim(phi.smooth.bs)[1], function(i){
  sort(phi.smooth.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))



df_phi = as.data.frame(cbind(Time = 1:dim(phi.ci.bs)[1], 
                             m.lower = phi.ci.bs[,1], 
                             m.upper = phi.ci.bs[,2], 
                             mu = rowMeans(phi.ci.bs),
                             true = true.para$phi.t.true[1:dim(phi.ci.bs)[1]]))

pdf("ci_phi.pdf")
s[[4]] <- ggplot(df_phi, aes(Time)) + 
  geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
              alpha=0.7)+
  geom_line(aes(y = mu), linetype = 'solid', size = 1,col='red') +
  geom_line(aes(y = true), linetype = 'solid', size = 1,col='blue') +
  labs(#title =expression(paste("'Pointwise Confidence Interval for ",phi,'(t)')),
       x = "Time", y = expression(paste(phi,'(t)'))) +theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))
dev.off()

grid.arrange(s[[1]],s[[2]],s[[4]],s[[3]],ncol=2,nrow=2)
## point-wise CI for del.t.lowess
del.lowess.bs <- est.bs$del.lowess.bs
#range(rhoH.bs)
del.t.ci.bs <- t(sapply(1:dim(del.lowess.bs)[1], function(i){
  sort(del.lowess.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

df_del = as.data.frame(cbind(Time = rep(1:(dim(del.t.ci.bs)[1]/10),10), 
                             m.lower = del.t.ci.bs [,1], 
                             m.upper = del.t.ci.bs [,2], 
                             mu = rowMeans(del.t.ci.bs ),
                             true = c(mortality_rate$mort_lowess_rate)))

pdf("ci_del_lowess.pdf")
p <- list()
for(i in 1:10){
  
  df_del_sub <- df_del[((i-1)*(dim(del.t.ci.bs)[1]/10)+1):(dim(del.t.ci.bs)[1]/10*i),]
  p[[i]] <- ggplot(df_del_sub, aes(Time)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu), linetype = 'solid', size = 1) +
    geom_line(aes(y = true), linetype = 'solid', size = 1,col='red') +
    labs(title = as.expression(bquote(paste("smoothed ",delta," for state"~.(i)))),
         x = "Time", y = expression(paste(delta,'(t)'))) +theme_bw()+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8))
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
             p[[10]],nrow=4,ncol=3)

dev.off()


##CI for del.t.lm

del.lm.bs <- est.bs$del.lm.bs

df_delta = as.data.frame(cbind(delta = del.lm.bs, delta.true = mortality_rate$mort_lm_rate,
                               del.mean = mean(del.lm.bs)))
pdf('hist_delta.pdf')
p <- list()
for(i in 1:10){
  
  df_delta_sub <- as.data.frame(cbind(delta = t(df_delta[i,]), delta.true = mortality_rate$mort_lm_rate[i],
                                      del.mean = mean(del.lm.bs)))
  colnames(df_delta_sub)[1] = 'delta'
  p[[i]] <- ggplot(df_delta_sub,aes(x=delta,y = ..density..))+
    geom_histogram(color="darkblue", 
                   fill="lightblue")+
    geom_density(color='darkred')+
    labs(title = as.expression(bquote(paste(bar(delta)," for state"~.(i)))),
         x = "Time", y = 'Density') +theme_bw()+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8)) +
    geom_vline(aes(xintercept = delta.true),col='darkgreen',size=1)+
    geom_vline(aes(xintercept = del.mean),col='red',size=1)
}
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
             p[[10]],nrow=4,ncol=3)

dev.off()

sd(del.lm.bs)
mean(del.lm.bs)
sort(del.lm.bs)[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]



## CI for delC
## under report one: (\theta+\gamma)At

theta.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$theta.hat.t)})

eta.bs <- sapply(1:sim.BS.N, function(i){return(theta.bs[,i]+est.bs$gamma.bs[i])}) 

delCt.bs <- eta.bs*At.bs

delCt.bs <- delCt.bs[-seq(94,940,94),] # remove the last day for each state
# to make the dimension is consistent for plot

delCt.ci.bs <- t(sapply(1:dim(delCt.bs)[1], function(i){
  sort(delCt.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))
delCt.true.smooth <- c(sapply(1:10,function(i){lowess(diff(dat.sim[[i]]$data$Ct[1:(dim(eta.bs)[1]/10)]),f=1/16)$y}))
#delCt.true <- c(sapply(1:10,function(i){diff(dat.sim[[i]]$data$Ct[1:(dim(eta.bs)[1]/10)])}))

df_delCt = as.data.frame(cbind(Time = rep(1:(dim(delCt.bs)[1]/10),10), m.lower = delCt.ci.bs[,1], 
                               m.upper = delCt.ci.bs[,2], mu = rowMeans(delCt.bs),
                               true = delCt.true.smooth))

pdf("ci_delCt.pdf")
p <- list()
for(i in 1:10){
  
  df_delCt_sub <- df_delCt[((i-1)*(dim(delCt.bs)[1]/10)+1):(dim(delCt.bs)[1]/10*i),]
  p[[i]] <- ggplot(df_delCt_sub, aes(Time)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu), linetype = 'solid', size = 1,col='red') +
    geom_line(aes(y = true), linetype = 'solid', size = 1,col='blue') +
    labs(#title = as.expression(bquote(paste('CI of ',Delta,'C(t) for state'~.(i)))),
         x = "Time", y = as.expression(bquote(paste(Delta,"C"[t],"-process"~.(i))))) +
    theme_bw()+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=18))
}

grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
             p[[10]],nrow=4,ncol=3)
dev.off()


