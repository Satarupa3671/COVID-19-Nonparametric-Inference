###Based on the output of the resi_boot.R, plots the CIs for different parameters for different states of the US. 

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
#source("./est_module_MN.R")
level <- 0.05
Nsim = 1000
load("../USA_data_processing/covid_df.Rda")

collect_parameters_boot = function(est.list, Nsim,level){
  ## collect gamma
  gamma.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$gamma.hat)
  ## collect rhoA
  rhoA.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$rhoA.hat)
  ## collect At
  At.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$A.hat.t)
  ## collect Ct
  Ct.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$cum.new.infect.t)
  ## point-wise CI for rhoH.t, phi.t
  rhoH.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$rhoH.final)
  
  ###theta.t
  theta.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$theta.hat.t)
  
  phi.smooth.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$est$phi.hat.smooth)
  
  ## point-wise CI for del.t.lowess
  del.lowess.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$mortality$mort_lowess_rate)
  
  ##CI for del.t.lm
  del.lm.bs <- sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_para$mortality$mort_lm_rate)
  ##CI for alpha.t
  ## collect alpha.t
  alpha.t.bs <- sapply(1:Nsim, function(i){
    (est.list[[i]]$boot_est$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:length(est.list[[i]]$boot_est$est_para$est$sqrt.alpha.kappa.t)])^2
  })
  
  sqrt.alpha.kappa.t.bs <- sapply(1:Nsim, function(i){
    est.list[[i]]$boot_est$est_para$est$sqrt.alpha.kappa.t
  })
  
  est.bs = list(del.lm.bs = del.lm.bs, del.lowess.bs  = del.lowess.bs , 
                sqrt.alpha.kappa.t.bs = sqrt.alpha.kappa.t.bs,
                phi.smooth.bs  = phi.smooth.bs , rhoH.bs = rhoH.bs, theta.bs = theta.bs,
                alpha.t.bs = alpha.t.bs,
                Ct.bs = Ct.bs ,At.bs= At.bs, gamma.bs = gamma.bs, rhoA.bs = rhoA.bs)
  save(est.bs,file = sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_est_bs.Rda", st, st, st))
  return(est.bs)
}

st_list = c("Arkansas", "Delaware", "Idaho", "Iowa", "Nebraska", "Ohio", "Oklahoma",
           "Pennsylvania", "South Dakota", "Tennessee", "Texas", "Utah", "Wisconsin")
for(st in st_list){
  data.st = data_st(st)
  load(sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_bootstrap_E.Rda", st, st, st))
  
  parameters_boot = collect_parameters_boot(est.list = est.list, Nsim = Nsim, level = level)  
  C_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$C)
  R.Reported_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$R.Reported)
  R_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$R)
  Q_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$Q)
  H_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$H)
  D_boot = sapply(1:Nsim, function(i) est.list[[i]]$boot_est$est_states$D)
  est_boot_data = list(C_boot = C_boot,
                       R.Reported_boot = R.Reported_boot ,
                       R_boot = R_boot,
                       Q_boot = Q_boot,
                       H_boot = H_boot,
                       D_boot = D_boot)
  
  save(est_boot_data,file = sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_est_boot_data.Rda", st, st, st))
}

for(st in st_list){ 
  load(sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_est_bs.Rda", st, st, st))
  load(sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_est_boot_data.Rda", st, st, st))
  load(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  
  data.st = data_st(st)
  
  n.t = est_st$est_para$est$n.t
  true.para <- list(gamma.true = est_st$est_para$est$gamma.hat,
                    rhoA.true = est_st$est_para$est$rhoA.hat,
                    rhoH.t.true = est_st$est_para$est$rhoH.final,
                    alpha.true = (est_st$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:(n.t-1)])^2,
                    delta.t.true = est_st$est_para$mortality$mort_lm_rate,
                    phi.t.true = est_st$est_para$est$phi.hat.smooth)
  
  
  Nsim = 1000
  level  = .05
  
  ## collect gamma
  
  File <- sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/gamma_ci.eps", st,  st, st)              
  
  if (file.exists(File)) stop(File, " already exists")
  dir.create(dirname(File), showWarnings = FALSE)
  
  df_gamma = as.data.frame(cbind(gamma = est.bs$gamma.bs, gamma.true = true.para$gamma.true,
                                 gamma.mean = mean(est.bs$gamma.bs)))
  p1 = ggplot(df_gamma, aes(x=gamma, y = ..density..))+
    geom_histogram(color="darkorange", fill= "darkorange1")+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',gamma)), 
      x = expression(gamma), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.y=element_blank(),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    geom_vline(aes(xintercept = gamma.true),col='blue',size=1)+
    geom_vline(aes(xintercept = gamma.mean),col='red',size=1.5)
  p1
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/gamma_ci.eps", st, st, st))
  
  
  ## collect rhoA
  df_rhoA = as.data.frame(cbind(rhoA = est.bs$rhoA.bs, rhoA.true =  true.para$rhoA.true,
                                rhoA.mean  = mean(est.bs$rhoA.bs)))
  p2 = ggplot(df_rhoA, aes(x=rhoA, y = ..density..))+
    geom_histogram(color="darkorange", fill= "darkorange1")+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',rho[A])), 
      x = expression(rho[A]), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.y=element_blank(),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26)) +
    geom_vline(aes(xintercept = rhoA.true),col='blue',size=1)+
    geom_vline(aes(xintercept = rhoA.mean),col='red',size=1.5)
  p2
 # { 
  #   rhoA = est.bs$rhoA.bs
  #   rhoA_updated = c(rep(unique(rhoA),length(rhoA)/4), 
  #                    rep(unique(rhoA)+0.0005,length(rhoA)/2),
  #                    rep(unique(rhoA)-0.0001,length(rhoA)/4)) 
  #   
  #   df_rhoA = as.data.frame(cbind(rhoA = rhoA_updated, rhoA.true =  true.para$rhoA.true,
  #                                 rhoA.mean  = mean(rhoA_updated)))
  #   p2 = ggplot(df_rhoA, aes(x=rhoA, y = ..density..))+
  #     geom_histogram(color="darkorange", fill= "darkorange1")+
  #     geom_density(color="darkred")+
  #     labs(#title = expression(paste('Sampling distribution of ',rho[A])), 
  #       x = expression(rho[A]), y = "Density") +theme_bw()+
  #     theme(axis.text=element_text(size=20),
  #           axis.text.y=element_blank(),
  #           axis.text.x = element_text(size=20,angle = 90, hjust = 1),
  #           axis.title=element_text(size=26)) +
  #     geom_vline(aes(xintercept = rhoA.true),col='blue',size=1)+
  #     geom_vline(aes(xintercept = rhoA.mean),col='red',size=1)
  #   p2
  # }  

  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/rhoA_ci.eps", st, st, st))
  
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
  
  
  p3 = ggplot(df_alpha, aes(x=alpha, y = ..density..))+
    geom_histogram(color="darkorange", fill= "darkorange1")+
    geom_density(color="darkred")+
    labs(#title = expression(paste('Sampling distribution of ',alpha)),
      x = expression(alpha), y = "Density") +theme_bw()+
    theme(axis.text=element_text(size=20),
          axis.text.y=element_blank(),
          axis.text.x = element_text(size=20,angle = 90, hjust = 1),
          axis.title=element_text(size=26)) +
    geom_vline(aes(xintercept = alpha.true),col='blue',size=1)+
    geom_vline(aes(xintercept = alpha.mean),col='firebrick2',size=1.5)
  p3
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/alpha_mean_ci.eps", st, st, st))
  # 
  
  
  ## point-wise CI for rhoH.t
  rhoH.bs <- est.bs$rhoH.bs
  # rhoH.bs = sapply(1:Nsim, function(i){
  #   n_eff = min(length(diff(est_boot_data$H_boot[,i])), length(est.bs$At.bs[,i]))
  #   pmin(.98,pmax(0,((mean(est.bs$gamma.bs)* (est.bs$At.bs[,i][1:n_eff] +
  #                                               est_boot_data$Q_boot[,i][1:n_eff]) -
  #                       diff(est_boot_data$H_boot[,i])[1:n_eff] -
  #                       est.bs$del.lowess.bs[,i][1:n_eff]*
  #                       est_boot_data$H_boot[,i][1:n_eff])/
  #                      est_boot_data$H_boot[,i][1:n_eff])))
  # })
  # rhoH.ci.bs <- t(sapply(1:dim(rhoH.bs)[1], function(i){
  #   sort(rhoH.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]}))


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
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+
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
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/rhoH_ci.eps", st,st, st))
  
  ## point-wise CI for theta.t
  theta.bs <- est.bs$theta.bs
  theta.ci.bs <- t(sapply(1:dim(theta.bs)[1], function(i){
    sort(theta.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]}))
  
  
  
  df_theta = as.data.frame(cbind(#Time = 1:dim(rhoH.ci.bs)[1], 
    m.lower = theta.ci.bs[,1], m.upper = theta.ci.bs[,2],
    mu = rowMeans(theta.ci.bs),
    theta.true = est_st$est_para$est$theta.hat.t[1:dim(theta.ci.bs)[1]]))
  df_theta$ddates = as.Date(data.st$dates[1:dim(theta.ci.bs)[1]], format = "%m/%d/%Y") 
  
  p5 = ggplot(df_theta, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+ 
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
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/theta_ci.eps", st, st, st))
  
  
  ## point-wise CI for del.t.lowess
  del.lowess.bs <- sapply(1:Nsim, function(i) est.bs$del.lowess.bs[,i])
  del.t.ci.bs <- abs(t(sapply(1:dim(del.lowess.bs)[1], function(i){
    sort(del.lowess.bs[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]})))
  
  df_del = as.data.frame(cbind(Time = 1:dim(del.t.ci.bs )[1], 
                               m.lower = del.t.ci.bs [,1], m.upper = del.t.ci.bs [,2], 
                               mu = rowMeans(del.t.ci.bs ) ,
                               del.true = est_st$est_para$mortality$mort_lowess_rate[ 1:dim(del.t.ci.bs)[1]]))
  df_del$ddates = as.Date(data.st$dates[1:dim(del.t.ci.bs)[1]], format = "%m/%d/%Y") 
  p6 = ggplot(df_del, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+
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
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/del_lowess_ci.eps", st,st, st))
  
  
  ###
  C <- est_boot_data$C_boot
  C.bs <- abs(t(sapply(1:dim(C)[1], function(i){
    sort(C[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]})))
  
  df_C = as.data.frame(cbind(Time = 1:dim(C.bs)[1], 
                             m.lower = C.bs[,1], m.upper = C.bs[,2], 
                             mu = rowMeans(C.bs) ,
                             C.true = est_st$est_states$C[1:dim(C.bs)[1]]))
  df_C$ddates = as.Date(data.st$dates[1:dim(C.bs )[1]], format = "%m/%d/%Y") 
  break.vec <- c(df_C$ddates[1],
                 seq(from = df_C$ddates[1] , to = df_C$ddates[length(df_C$ddates)], 
                     by="month"), df_C$ddates[length(df_C$ddates)])
  
  p7 = ggplot(df_C, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+
    geom_line(aes(y = mu), col = 'red', linetype = 'solid', size = 1) +
    geom_line(aes(y = C.true), linetype = 'solid', col ='blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",delta,'(t)')),
      x = "Time", y = expression(paste(C,'(t)'))) +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26 ))
  
  p7
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/C.eps", st, st, st))
  
  
  ####
  R <- est_boot_data$R_boot
  R.bs <- abs(t(sapply(1:dim(R)[1], function(i){
    sort(R[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]})))
  
  df_R = as.data.frame(cbind(Time = 1:dim(R.bs)[1], 
                             m.lower = R.bs[,1], m.upper = R.bs[,2], 
                             mu = rowMeans(R.bs) ,
                             R.true = est_st$est_states$R[1:dim(R.bs)[1]]))
  df_R$ddates = as.Date(data.st$dates[1:dim(R.bs )[1]], format = "%m/%d/%Y") 
  p8 = ggplot(df_R, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+
    geom_line(aes(y = mu), col = 'red', linetype = 'solid', size = 1) +
    geom_line(aes(y = R.true), linetype = 'solid', col ='blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",delta,'(t)')),
      x = "Time", y = expression(paste(R,'(t)'))) +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26 ))
  
  p8
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/R.eps", st, st, st))
  
  
  ####
  R <- est_boot_data$R.Reported_boot
  R_repo.bs <- abs(t(sapply(1:dim(R)[1], function(i){
    sort(R[i,])[round(c((Nsim+1)*level/2, (Nsim+1)*(1-level/2)))]})))
  
  df_R_repo = as.data.frame(cbind(Time = 1:dim(R_repo.bs)[1], 
                                  m.lower = R_repo.bs[,1], m.upper = R_repo.bs[,2], 
                                  mu = rowMeans(R_repo.bs) ,
                                  R_repo.true = est_st$est_states$R.Reported[1:dim(R.bs)[1]]))
  df_R_repo$ddates = as.Date(data.st$dates[1:dim(R_repo.bs )[1]], format = "%m/%d/%Y") 
  p9 = ggplot(df_R_repo, aes(ddates)) + 
    geom_ribbon(aes(ymin=m.lower,ymax=m.upper), fill = "darkorange1")+
    geom_line(aes(y = mu), col = 'red', linetype = 'solid', size = 1) +
    geom_line(aes(y = R_repo.true), linetype = 'solid', col ='blue', size = 1) +
    labs(#title =expression(paste("'Pointwise Confidence Interval for ",delta,'(t)')),
      x = "Time", y = expression(paste(R,'(t)'))) +theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26 ))
  
  p9
  ggsave(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/R_repo.eps", st, st, st))
  
  
  
  
  sink(sprintf("../Results for USA States/%s/%s_updated_boot_final/CI_%s/CIs.txt", st,st, st))
  
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
