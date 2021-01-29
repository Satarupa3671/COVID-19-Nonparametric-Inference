
##Particularly for the state Minnesota (MN) --to avoid the recording discrepancy in the observed data 
### The function estimation_module -- takes in any observed data, NumBlock and wd.length and returns the estimated parameters, as well as the estimates for the states C,Q,R,D,H. Based on the input argument boot == FALSE (based on the observed sample, not the bootstrapped sample), it saves the ggplot of the profile loss vs gamma.grid and profile loss vs rhoA.grid.
###The function plot_for_estimation --  runs the function estimation_module for boot == FALSE (observed data only) and saves the estimated parameters etc.


rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
#library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(scales)

source("./fitting_functions.R")
source("./est_plot_func.R")
load("../USA_data_processing/covid_df.Rda")

##data_st is the function to read the data for a particular state "st" and return the smoothed trajectories of (C,D,H,Q,R,DelC,DelH,DelR,DelQ,DelD,Test,Kappa) as a list
data_st = function(st){
  covid.state = subset(covid_final, state == st)
  ##which columns contain NA values?
  colSums(covid.state[,c(3:6,9:ncol(covid.state))]) ##hospitalizedCurrently , hospitali+   ## from which date  is the data available 
  #date_notNA1 = max(max(covid.state$date[is.na(covid.state$recovered)]),
  #                max(covid.state$date[is.na(covid.state$hospitalizedCurrently)]))
  date_notNA = max(covid.state$date[is.na(covid.state$recovered)])
  
  date_notNA1 = covid.state$date[is.na(covid.state$hospitalizedCurrently)]
  covid.state = subset(covid.state, date > date_notNA)
  covid.state = subset(covid.state, !(covid.state$date  %in% date_notNA1) )
  delC = lowess(diff(covid.state$total_cases), f = bw_lowess)$y
  delR = lowess(diff(covid.state$recovered), f = bw_lowess)$y
  delH = lowess(diff(covid.state$hospitalizedCurrently), f = bw_lowess)$y
  delD = lowess(diff(covid.state$total_deaths), f = bw_lowess)$y
  
  Q = lowess((covid.state$total_cases - covid.state$recovered - covid.state$total_deaths - covid.state$hospitalizedCurrently), f = bw_lowess)$y
  H = lowess(covid.state$hospitalizedCurrently, f = bw_lowess)$y
  delTests = lowess((diff(covid.state$positive + covid.state$negative)), f = bw_lowess)$y
  TtoH_deno = lowess(covid.state$hospitalizedCurrently[-length(covid.state$hospitalizedCurrently)], f = bw_lowess)$y
  TtoH = delTests/TtoH_deno
  kappa = (100 + 100*covid.state$average)/100
  pop19 = covid.state$pop19
  var = list(delC =  delC, 
             delR = delR, delH = delH, delD = delD,
             Q = Q, H = H, TtoH = TtoH, pop19= pop19,
             R = lowess(covid.state$recovered, f = bw_lowess)$y, 
             D = lowess(covid.state$total_deaths, f = bw_lowess)$y,
             C = lowess(covid.state$total_cases, f = bw_lowess)$y,
             Test = lowess((covid.state$positive + covid.state$negative), f= bw_lowess)$y,
             kappa = kappa)
  
  raw_data = list(delC = diff(covid.state$total_cases) , 
                  delR = diff(covid.state$recovered), 
                  delH = diff(covid.state$hospitalizedCurrently), 
                  delD = diff(covid.state$total_deaths),
                  Q = (covid.state$total_cases - covid.state$recovered - 
                         covid.state$total_deaths - covid.state$hospitalizedCurrently),
                  H = covid.state$hospitalizedCurrently, 
                  TtoH = (diff(covid.state$positive + covid.state$negative))/covid.state$hospitalizedCurrently[1:length(diff(covid.state$positive + covid.state$negative))],
                  R = covid.state$recovered, D = covid.state$total_deaths,
                  C = covid.state$total_cases,
                  Test = (covid.state$positive + covid.state$negative))
  n = length(var$Q)
  return(list(st = st,var = var, raw_data = raw_data, n=n, dates = covid.state$date))
}

estimated_states = function(dat.int, true.para,Time){
  
  A <- dat.int$At
  S <- c(dat.int$St,rep(0,Time-2))
  R <- c(dat.int$Rt,rep(0,Time-2))
  R.Reported <- c(dat.int$Rt,rep(0,Time-2))
  Q <- c(dat.int$Qt,rep(0,Time-2))
  H <- c(dat.int$Ht,rep(0,Time-2))
  D <- c(dat.int$Dt,rep(0,Time-2))
  C <- c(dat.int$Ct,rep(0,Time-2))
  Test <- dat.int$T.t
  sqrt.alpha.kappa.t <- dat.int$sqrt.alpha.kappa.t
  
  gamma <- true.para$gamma.true
  rhoA <- true.para$rhoA.true
  rhoH.t <- true.para$rhoH.t.true
  alpha <- true.para$alpha.true
  delta.t <- true.para$delta.t.true
  phi.t <- true.para$phi.t.true
  
  theta.t = dat.int$theta.hat.t
  for (t in 1:(Time-2)){
    DelC.t = (theta.t[t] + gamma)*A[t]
    #C[t+1] = pmin((C[t]+ DelC.t),Max.pop)
    C[t+1] = abs(C[t]+ DelC.t)
    
    DelQ.t = theta.t[t]*A[t] - (gamma + rhoA)*Q[t]
    #Q[t+1] = pmin((Q[t] + DelQ.t), Max.pop)
    Q[t+1] = abs(Q[t] + DelQ.t)
    
    DelH.t = gamma*A[t] + gamma*Q[t] - rhoH.t[t]*H[t] - delta.t[t]*H[t]
    #H[t+1] = pmin((H[t] + DelH.t), C[t])
    H[t+1] = abs(H[t] + DelH.t)
    
    DelRH.t = rhoH.t[t]*H[t] 
    DelRQ.t = rhoA*Q[t]
    DelRA.t = rhoA*A[t]
    
    DelR.t = DelRH.t + DelRQ.t +DelRA.t
    #R[t+1] = pmin((R[t] + DelR.t),C[t])
    R[t+1] = abs(R[t] + DelR.t)
    
    DelR.Reported.t = DelRH.t + DelRQ.t #+DelRA.t
    #R[t+1] = pmin((R[t] + DelR.t),C[t])
    R.Reported[t+1] = abs(R.Reported[t] + DelR.Reported.t)
    
    DelD.t = delta.t[t]*H[t]
    #D[t+1] = pmin((D[t] + DelD.t),C[t])
    D[t+1] = abs(D[t] + DelD.t)
    
    DelS.t =  -(sqrt.alpha.kappa.t[t])^2*A[t] * S[t]/(S[t]+R[t]+A[t])
    #S[t+1] = pmin((S[t] + DelS.t),Max.pop)
    S[t+1] = abs(S[t] + DelS.t)
  }
  
  # C = lowess(C, f = bw_lowess)$y
  # R = lowess(R, f = bw_lowess)$y
  # Q = lowess(Q, f = bw_lowess)$y
  # H = lowess(H, f = bw_lowess)$y
  # D = lowess(D, f = bw_lowess)$y
  return(list(A = A, S = S, R = R, Q = Q, H = H, D = D, C = C,
              R.Reported = R.Reported, Test = Test,
              sqrt.alpha.kappa.t = sqrt.alpha.kappa.t,
              theta.t = theta.t))
}


estimation_module = function(dat, numBlock, wd.length, boot = TRUE){
  
  n = dat$n
  shift =  ((n-1) - wd.length)/(numBlock - 1) #10
  TimeBlock = NULL
  for(j in 1:numBlock){
    TimeBlock[[j]] = (1 + shift * (j-1)) : min(shift * (j-1) + wd.length,(n-1))
  }
  
  gamma.grid = seq(0.0001,0.03,0.0002)
  rhoA.grid = seq(0.002,0.2,0.002)
  eval.mat = matrix(0,length(gamma.grid),length(rhoA.grid))
  
  for(i in 1:length(gamma.grid)){
    for(j in 1:length(rhoA.grid)){
      eval.mat[i,j] = objective.grA(dat$var,c(gamma.grid[i],rhoA.grid[j]),TimeBlock)
    }
  }
  
  eval.gamma = apply(eval.mat,1,min)
  eval.rhoA = apply(eval.mat,2,min)
  
  g.index = which(eval.gamma == min(eval.gamma )) # which.min(eval.gamma)
  rA.index = which(eval.rhoA == min(eval.rhoA )) # which.min(eval.rhoA)
  
  gamma.hat = gamma.grid[g.index][1]
  rhoA.hat = rhoA.grid[rA.index]
  
  if(boot == FALSE){
    
    df_gamma = as.data.frame(cbind(
      gamma.grid = gamma.grid, eval.gamma = eval.gamma))
    
    ##ggplot of profile loss as a func of gamma 
    File <- sprintf("../Results for USA States/%s/%s_plot_gamma.pdf", dat$st, dat$st)              
    #if (file.exists(File)) stop(File, " already exists")
    dir.create(dirname(File), showWarnings = FALSE)
    
    ggplot() + geom_line(data = df_gamma,aes(x = gamma.grid,y = eval.gamma)) +
      xlab(expression(gamma)) +
      ylab("Profile Loss")+
      theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14))+
      theme(legend.position= "top",
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "transparent"), 
            legend.text=element_text(size=14))
    ggsave(sprintf("../Results for USA States/%s/%s_plot_gamma.pdf",  dat$st,  dat$st))
    
    
    ##ggplot of profile loss as a func of rhoA
    df_rhoA = as.data.frame(cbind(rhoA.grid = rhoA.grid, eval.rhoA = eval.rhoA))
    ggplot() + geom_line(data = df_rhoA,aes(x = rhoA.grid,y = eval.rhoA)) +
      xlab(expression(rho[A])) +
      ylab("Profile Loss")+
      theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14))+
      theme(legend.position= "top",
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "transparent"), 
            legend.text=element_text(size=14))
    
    ggsave(sprintf("../Results for USA States/%s/%s_plot_rhoA.pdf",  dat$st,  dat$st))
    
  }
  
  mort_lowess_rate = lowess((dat$var$delD/dat$var$H[-length(dat$var$H)]),f=bw_lowess)$y
  mort_lm = lm(dat$var$delD ~ dat$var$H[-length(dat$var$H)] -1)
  mort_lm_rate = mort_lm$coefficients[1]
  
  est1 = estimate_parameters(gamma.hat = gamma.hat, rhoA.hat = rhoA.hat, 
                             var = dat$var, shift = shift, 
                             numBlock = numBlock, 
                             wd.length = wd.length, TimeBlock = TimeBlock)
  
  if(boot == FALSE){  true.para = list(gamma.true = est1$gamma.hat,
                                       rhoA.true = est1$rhoA.hat,
                                       rhoH.t.true = est1$rhoH.final,
                                       alpha.true = (mean(est1$sqrt.alpha.kappa.t/dat$kappa[1:length(est1$sqrt.alpha.kappa.t)]))^2,
                                       delta.t.true = mort_lowess_rate,
                                       phi.t.true = est1$phi.hat.smooth)
  }else {true.para = true.para.st}
  
  dat.int <- list(At = est1$A.hat.t, 
                  sqrt.alpha.kappa.t = est1$sqrt.alpha.kappa.t,
                  theta.hat.t = est1$theta.hat.t,
                  Rt = dat$var$R[1],
                  Qt = dat$var$C[1] - dat$var$R[1] - dat$var$D[1] - dat$var$H[1],
                  Ht = dat$var$H[1],
                  Dt = dat$var$D[1],
                  Ct = dat$var$C[1],
                  T.t = dat$var$Test,
                  kappa = dat$var$kappa)
  dat.int$St = unique(dat$var$pop19) -  dat.int$At[1] -  dat.int$Rt - dat.int$Qt - dat.int$Ht - dat.int$Dt - dat.int$Ct
  
  est_states = estimated_states(dat.int = dat.int, true.para = true.para, Time = n)
  
  return(list(est_para = list(est = est1, 
                              mortality = list(mort_lm_rate = mort_lm_rate, mort_lowess_rate = mort_lowess_rate)),
              est_states = est_states, 
              tuning = list(numBlock = numBlock, wd.length = wd.length, TimeBlock = TimeBlock)
              
  ))
}

####

plot_for_estimation = function(st, numBlock, wd.length){
  data.st = data_st(st)
  est_st = estimation_module(dat = data.st, numBlock = numBlock , wd.length = wd.length, boot = FALSE)
  #load(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  
  sink(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  save(est_st, file = sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  sink()
  plot_est(var = data.st$var, raw_data = data.st$raw_data, est = est_st, dates = data.st$dates)
  
  # numBlock = est_st$tuning$numBlock
  # TimeBlock = est_st$tuning$TimeBlock
  # pdf(sprintf("./%s/%s_plot_diagnostics.pdf", st, st))
  # for(j in 1:numBlock){
  #   diagnostics(whichBlock = j, TimeBlock = TimeBlock , est = est1 , 
  #               var = data.st$var)
  # }
  # dev.off()
}
# 
# st = "Minnesota"
# numBlock = 12 #opt_wd$numBl_st[opt_wd$states == st]
# wd.length = 29 #opt_wd$wdL_st[opt_wd$states == st]
# plot_for_estimation(st, numBlock, wd.length)

