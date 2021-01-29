###Based on the output of the resi_boot.R, plots the CIs for different parameters across different states. 
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(scales)
source("./est_module.R")

st_list = c("Arkansas", "Delaware", "Idaho", "Iowa", 
            "Nebraska","Ohio", "Oklahoma", "Pennsylvania", "South Dakota" ,
            "Tennessee","Texas", "Wisconsin")

load("../USA_data_processing/covid_df.Rda")

for(st in st_list){ 
  load(sprintf("../Results for USA States/%s/%s_est_bs_1120.Rda", st, st))
  load(sprintf("../Results for USA States/%s/%s_est_boot_data_1120.Rda", st, st))
  load(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  
  # load(sprintf("../nonparametric_est_data_analysis_Dec15/Utah/Utah_est_bs_1120.Rda"))
  # load(sprintf("../nonparametric_est_data_analysis_Dec15/Utah/Utah_est_boot_data_1120.Rda"))
  # load(sprintf("../nonparametric_est_data_analysis_Dec15/Utah/Utah_est_st.Rda"))
  # 
  data.st = data_st(st)
  
  n.t = est_st$est_para$est$n.t
  true.para <- list(gamma.true = est_st$est_para$est$gamma.hat,
                    rhoA.true = est_st$est_para$est$rhoA.hat,
                    rhoH.t.true = est_st$est_para$est$rhoH.final,
                    alpha.true = (est_st$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:(n.t-1)])^2,
                    delta.t.true = est_st$est_para$mortality$mort_lm_rate,
                    phi.t.true = est_st$est_para$est$phi.hat.smooth)
  
  
  Nsim = 1050
  level  = .05
  ## collect gamma
  
  File <- sprintf("../Results for USA States/%s/CI_%s/gamma_ci.pdf", st, st)              
  
  #if (file.exists(File)) stop(File, " already exists")
  dir.create(dirname(File), showWarnings = FALSE)
  
  df_gamma = as.data.frame(cbind(gamma = est.bs$gamma.bs, gamma.true = true.para$gamma.true,
                                 gamma.mean = mean(est.bs$gamma.bs)))
  p1 = ggplot(df_gamma, aes(x=gamma, y = ..density..))+
    geom_histogram(color="darkorange", fill= "darkorange1", alpha=.7) + # "lightblue")+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',gamma)), 
      x = expression(gamma), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    geom_vline(aes(xintercept = gamma.true),col='blue',size=1)+
    geom_vline(aes(xintercept = gamma.mean),col='red',size=.6)
  p1
  ggsave(sprintf("../Results for USA States/%s/CI_%s/gamma_ci.pdf", st, st))
  
  
  ## collect rhoA
  df_rhoA = as.data.frame(cbind(rhoA = est.bs$rhoA.bs, rhoA.true =  true.para$rhoA.true,
                                rhoA.mean  = mean(est.bs$rhoA.bs)))
  p2 = ggplot(df_rhoA, aes(x=rhoA, y = ..density..))+
    geom_histogram(color="darkorange", fill= "darkorange1", alpha=.7)+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',rho[A])), 
      x = expression(rho[A]), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26)) +
    geom_vline(aes(xintercept = rhoA.true),col='blue',size=1)+
    geom_vline(aes(xintercept = rhoA.mean),col='red',size=.4)
  p2
  ggsave(sprintf("../Results for USA States/%s/CI_%s/rhoA_ci.pdf", st, st))
  
  ####pointwise ci for alpha.t
  
  ## collect alpha
  
  
  alpha.t.bs <- est.bs$alpha.t.bs
  alpha.ci.bs <- t(sapply(1:dim(alpha.t.bs)[1], function(i){
    sort(alpha.t.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]}))
  
  df_alpha.t = as.data.frame(cbind(#Time = 1:dim(alpha.ci.bs)[1],
    m.lower = alpha.ci.bs[,1], m.upper = alpha.ci.bs[,2], 
    mu = rowMeans(alpha.ci.bs),
    alpha.true = true.para$alpha.true[ 1:dim(alpha.ci.bs)[1]]))
  df_alpha.t$ddates = as.Date(data.st$dates[1:dim(alpha.ci.bs)[1]], format = "%m/%d/%Y") 
  # ggplot(df_alpha.t, aes(ddates)) + 
  #   geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
  #               alpha=0.7)+
  #   geom_line(aes(y = mu), linetype = 'solid', size = 1) +
  #   geom_line(aes(y = alpha.true), linetype = 'solid', col = 'red', size = 1) +
  #   labs(#title =expression(paste("'Pointwise Confidence Interval for ",alpha,'(t)')),
  #     x = "Time", y = expression((paste(alpha,'(t)')))) +theme_bw()+
  #   scale_x_date(breaks=date_breaks("1 months"),
  #                labels = date_format("%b %y")) +
  #   theme(axis.text=element_text(size=12), 
  #         axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
  #         axis.title=element_text(size=12))
  # 
  # ggsave(sprintf("./%s/CI_%s/alpha_ci.pdf", st, st))
  
  
  ###alpha.mean
  alpha_est_boot = matrix(0,dim(est.bs$alpha.t.bs)[1],Nsim)
  for(rep in 1:Nsim){
    for(time in 1:dim(est.bs$alpha.t.bs)[1]){
      alpha_est_boot[time, rep] = est.bs$alpha.t.bs[time, rep]*data.st$var$kappa[time]/mean(data.st$var$kappa)
    }
  }
  
  
  df_alpha = as.data.frame(cbind( alpha = colMeans(alpha_est_boot),
                                  alpha.true = mean(true.para$alpha.true),
                                  alpha.mean = mean(alpha_est_boot)))
  
  
  # df_alpha = as.data.frame(cbind(alpha = rowMeans(est.bs$alpha.t.bs), alpha.true = mean(true.para$alpha.true),
  #                                alpha.mean  = mean(est.bs$alpha.t.bs)))
  # 
  # 
  p3 = ggplot(df_alpha, aes(x=alpha, y = ..density..))+
   # geom_histogram(color="darkblue", fill="lightblue")+
    geom_histogram(color="darkorange", fill= "darkorange1", alpha=.7)+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',alpha)),
      x = expression(alpha), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26)) +
    geom_vline(aes(xintercept = alpha.true),col='blue',size=1)+
    geom_vline(aes(xintercept = alpha.mean),col='firebrick2',size=1)
  p3
  ggsave(sprintf("../Results for USA States/%s/CI_%s/alpha_mean_ci.pdf", st, st))
  # 
  
  
   ## point-wise CI for rhoH.t
  rhoH.bs <- est.bs$rhoH.bs
  #range(rhoH.bs)
  rhoH.ci.bs <- t(sapply(1:dim(rhoH.bs)[1], function(i){
    sort(rhoH.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]}))
  

  
  
  df_rhoH = as.data.frame(cbind(#Time = 1:dim(rhoH.ci.bs)[1], 
    m.lower = rhoH.ci.bs[,1], m.upper = rhoH.ci.bs[,2],
    mu = rowMeans(rhoH.ci.bs),
    rhoH.true = est_st$est_para$est$rhoH.final[1:dim(rhoH.ci.bs)[1]]))
  df_rhoH$ddates = as.Date(data.st$dates[1:dim(rhoH.ci.bs)[1]], format = "%m/%d/%Y") 
  
  break.vec <- c(df_rhoH$ddates[1],
                 seq(from = df_rhoH$ddates[1] , to = df_rhoH$ddates[length(df_rhoH$ddates)], 
                     by="month"), df_rhoH$ddates[length(df_rhoH$ddates)])
  
  p4 = ggplot(df_rhoH, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu),col ='red', linetype = 'solid', size = 1) +
    geom_line(aes(y = rhoH.true), linetype = 'solid', col = 'blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",rho[H],'(t)')),
      x = "Time", y = expression(paste(rho[H],'(t)')))  +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))
  
  p4
  ggsave(sprintf("../Results for USA States/%s/CI_%s/rhoH_ci.pdf", st, st))
  
  
  
  ## point-wise CI for theta.t
  theta.bs <- est.bs$theta.bs
  #range(rhoH.bs)
  theta.ci.bs <- t(sapply(1:dim(theta.bs)[1], function(i){
    sort(theta.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]}))
  
  
  
  df_theta = as.data.frame(cbind(#Time = 1:dim(rhoH.ci.bs)[1], 
    m.lower = theta.ci.bs[,1], m.upper = theta.ci.bs[,2],
    mu = rowMeans(theta.ci.bs),
    theta.true = est_st$est_para$est$theta.hat.t[1:dim(theta.ci.bs)[1]]))
  df_theta$ddates = as.Date(data.st$dates[1:dim(theta.ci.bs)[1]], format = "%m/%d/%Y") 
  
  p5 = ggplot(df_theta, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=0.7)+
    geom_line(aes(y = mu),col ='red', linetype = 'solid', size = 1) +
    geom_line(aes(y = theta.true), linetype = 'solid', col = 'blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",rho[H],'(t)')),
      x = "Time", y = expression(paste(theta,'(t)')))  +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))
  
  p5
  ggsave(sprintf("../Results for USA States/%s/CI_%s/theta_ci.pdf", st, st))
  
  
  ## point-wise CI for del.t.lowess
  del.lowess.bs <- est.bs$del.lowess.bs
  #range(rhoH.bs)
  del.t.ci.bs <- abs(t(sapply(1:dim(del.lowess.bs)[1], function(i){
    sort(del.lowess.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]})))
  
  df_del = as.data.frame(cbind(Time = 1:dim(del.t.ci.bs )[1], 
                               m.lower = del.t.ci.bs [,1], m.upper = del.t.ci.bs [,2], 
                               mu = rowMeans(del.t.ci.bs ) ,
                               del.true = est_st$est_para$mortality$mort_lowess_rate[ 1:dim(del.t.ci.bs)[1]]))
  df_del$ddates = as.Date(data.st$dates[1:dim(del.t.ci.bs)[1]], format = "%m/%d/%Y") 
  p6 = ggplot(df_del, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1",
                alpha=.7)+
    geom_line(aes(y = mu), col = 'red', linetype = 'solid', size = 1) +
    geom_line(aes(y = del.true), linetype = 'solid', col ='blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",delta,'(t)')),
      x = "Time", y = expression(paste(delta,'(t)'))) +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26 ))
  
  p6 
  ggsave(sprintf("../Results for USA States/%s/CI_%s/del_lowess_ci.pdf", st, st))
  
  sink(sprintf("../Results for USA States/%s/CI_%s/CIs.txt", st, st))
  
  print(paste("ci for gamma = ",sort(est.bs$gamma.bs)[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]))
  print(paste("mean of gamma.hat=", mean(est.bs$gamma.bs)))
  print(paste("sd of gamma.hat = ", sd(est.bs$gamma.bs)))
  
  print(paste("ci for rhoA = ",sort(est.bs$rhoA.bs)[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]))
  print(paste("mean of rhoA.hat=", mean(est.bs$rhoA.bs)))
  print(paste("sd of rhoA.hat = ", sd(est.bs$rhoA.bs)))
  
  print(paste("ci for alpha.mean = ",
              sort(df_alpha$alpha)[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]))
  print(paste("mean of alpha.mean=", mean(df_alpha$alpha)))
  print(paste("sd of alpha.mean = ", sd(df_alpha$alpha)))
  sink()
  
  
}
  
  