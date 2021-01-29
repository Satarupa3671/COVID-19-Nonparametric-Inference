###generates the bootstrapped samples based on the observed and estimated states C,Q,R,H,D.
##calls the estimation module function for the argument boot == TRUE. 
###Saves the estimated parameters and the estimated states C,Q,R,H,D (As an output of the function estimation_module) for each bootstrap sample into .Rda files.
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(parallel)

source("./fitting_functions.R")
source("./est_plot_func.R")
source("./est_module.R")

load("../USA_data_processing/covid_df.Rda")
load("./Utah/Utah_est_st.Rda")

resi_bootstrap = function(observed, estimated, Time){
  
  Max.pop = unique(observed$pop19)
  resi.Ct = observed$C[1:Time-1] - estimated$C
  resi.Rt = observed$R[1:Time-1] - estimated$R
  resi.Qt = observed$Q[1:Time-1] - estimated$Q
  resi.Ht = observed$H[1:Time-1] - estimated$H
  resi.Dt = observed$D[1:Time-1] - estimated$D
  
  resi.Ct.smooth = lowess(resi.Ct,f = 1/16)$y
  resi.Ct.adj = resi.Ct - resi.Ct.smooth
  resi.Ct.adj.samp = sample(resi.Ct.adj, replace = TRUE)
  Ct.adj.boot = abs(resi.Ct.smooth + resi.Ct.adj.samp  + estimated$C)
  
  
  resi.Rt.smooth = lowess(resi.Rt,f = 1/16)$y
  resi.Rt.adj = resi.Rt - resi.Rt.smooth
  resi.Rt.adj.samp = sample(resi.Rt.adj, replace = TRUE)
  Rt.adj.boot = abs(resi.Rt.smooth + resi.Rt.adj.samp  + estimated$R)
  
  
  resi.Qt.smooth = lowess(resi.Qt,f = 1/16)$y
  resi.Qt.adj = resi.Qt - resi.Qt.smooth
  resi.Qt.adj.samp = sample(resi.Qt.adj, replace = TRUE)
  Qt.adj.boot = abs(resi.Qt.smooth + resi.Qt.adj.samp  + estimated$Q)
  
  
  resi.Ht.smooth = lowess(resi.Ht,f = 1/16)$y
  resi.Ht.adj = resi.Ht - resi.Ht.smooth
  resi.Ht.adj.samp = sample(resi.Ht.adj, replace = TRUE)
  Ht.adj.boot = abs(resi.Ht.smooth + resi.Ht.adj.samp  + estimated$H)
  
  
  resi.Dt.smooth = lowess(resi.Dt,f = 1/16)$y
  resi.Dt.adj = resi.Dt - resi.Dt.smooth
  resi.Dt.adj.samp = sample(resi.Dt.adj, replace = TRUE)
  Dt.adj.boot = abs(resi.Dt.smooth + resi.Dt.adj.samp  + estimated$D)
  
  A = estimated$A
  St.adj.boot = abs(rep(Max.pop,length(A)) - A -  Rt.adj.boot[1:(Time-1)] - Qt.adj.boot[1:(Time-1)] - 
                      Ht.adj.boot[1:(Time-1)] - Dt.adj.boot[1:(Time-1)] - Ct.adj.boot[1:(Time-1)])
  
  delC = lowess(diff(Ct.adj.boot), f = 1/16)$y
  delR = lowess(diff(Rt.adj.boot), f= 1/16)$y
  delH = lowess(diff(Ht.adj.boot), f = 1/16)$y
  delD = lowess(diff(Dt.adj.boot), f = 1/16)$y
  delTests = lowess((diff(observed$Test)), f = 1/16)$y
  TtoH_deno = lowess((Ht.adj.boot), f = 1/16)$y
  TtoH = delTests/TtoH_deno[1:length(delTests)]
  
  var = list(delC = delC, delR = delR, delH = delH, delD = delD, 
             A = A,
             Q = Qt.adj.boot[1:(Time -1)], H = Ht.adj.boot[1:(Time -1)], TtoH = TtoH[1:(Time -1)],
             S = St.adj.boot[1:(Time -1)], R = Rt.adj.boot[1:(Time -1)], D = Dt.adj.boot[1:(Time -1)],
             C = Ct.adj.boot[1:(Time -1)], Test = observed$Test[1:(Time -1)],
             pop19 = observed$pop19,
             kappa= observed$kappa)
  
 # dat.sim <- data.frame(At=A) %>%
 #   mutate(St = St.adj.boot , Rt = Rt.adj.boot[1:(Time -1)],
  #         Qt = Qt.adj.boot[1:(Time -1)] ,Ht = Ht.adj.boot[1:(Time -1)],
   #        Dt = Dt.adj.boot[1:(Time -1)] ,Ct = Ct.adj.boot[1:(Time -1)],
    #       Tt = observed$Test[1:(Time -1)])
  
  return(list(var = var, theta = estimated$theta.t, n = Time))
}


ncore <-5
Nsim = 5
time <- proc.time()



###read all the variables
st = "Utah"
data.st = data_st(st)
observed_states = data.st
estimated = est_st$est_states

true.para.st <- list(gamma.true = est_st$est_para$est$gamma.hat,
                  rhoA.true = est_st$est_para$est$rhoA.hat,
                  rhoH.t.true = est_st$est_para$est$rhoH.final,
                  alpha.true = (mean(est_st$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:length(est_st$est_para$est$sqrt.alpha.kappa.t)]))^2,
                  delta.t.true = est_st$est_para$mortality$mort_lowess_rate,
                  phi.t.true = est_st$est_para$est$phi.hat.smooth)

est.list <- mclapply(1:Nsim, function(i){
  boot_sample = resi_bootstrap(observed = observed_states$var,
                               estimated =  estimated, Time = data.st$n)
  

  boot_est = estimation_module(dat = boot_sample, 
                           numBlock = est_st$tuning$numBlock , 
                           wd.length = est_st$tuning$wd.length,
                            boot = TRUE) 
  return(boot_est = boot_est)
}, mc.cores = ncore)




level <- 0.05

collect_parameters_boot = function(est.list, Nsim,level){
  ## collect gamma
  gamma.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$gamma.hat)
  ## collect rhoA
  rhoA.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$rhoA.hat)
  ## collect At
  At.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$A.hat.t)
  ## collect Ct
  Ct.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$cum.new.infect.t)
  ## point-wise CI for rhoH.t, phi.t
  rhoH.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$rhoH.final)
  
  ###theta.t
  theta.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$theta.hat.t)
  
  phi.smooth.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$est$phi.hat.smooth)
  
  ## point-wise CI for del.t.lowess
  del.lowess.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$mortality$mort_lowess_rate)
  
  ##CI for del.t.lm
  del.lm.bs <- sapply(1:Nsim, function(i) est.list[[i]]$est_para$mortality$mort_lm_rate)
  ##CI for alpha.t
  ## collect alpha.t
  alpha.t.bs <- sapply(1:Nsim, function(i){
    (est.list[[i]]$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:length(est.list[[i]]$est_para$est$sqrt.alpha.kappa.t)])^2
  })
  
  est.bs = list(del.lm.bs = del.lm.bs, del.lowess.bs  = del.lowess.bs , 
                phi.smooth.bs  = phi.smooth.bs , rhoH.bs = rhoH.bs, theta.bs = theta.bs,
                alpha.t.bs = alpha.t.bs,
                Ct.bs = Ct.bs ,At.bs= At.bs, gamma.bs = gamma.bs, rhoA.bs = rhoA.bs)
  save(est.bs,file = "Utah_est_bs_1120.Rda")
  return(est.bs)
}

parameters_boot = collect_parameters_boot(est.list = est.list, Nsim = Nsim, level = level)  


C_boot = sapply(1:Nsim, function(i) est.list[[i]]$est_states$C)
R_boot = sapply(1:Nsim, function(i) est.list[[i]]$est_states$R)
Q_boot = sapply(1:Nsim, function(i) est.list[[i]]$est_states$Q)
H_boot = sapply(1:Nsim, function(i) est.list[[i]]$est_states$H)
D_boot = sapply(1:Nsim, function(i) est.list[[i]]$est_states$D)
est_boot_data = list(C_boot = C_boot,
                     R_boot = R_boot,
                     Q_boot = Q_boot,
                     H_boot = H_boot,
                     D_boot = D_boot)

save(est_boot_data,file = "Utah_est_boot_data_1120.Rda")




