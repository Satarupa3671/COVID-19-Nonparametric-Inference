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
#source("./est_module_MN.R")
load("../USA_data_processing/covid_df.Rda")
bw_lowess = 1/16

##This function jointly resamples from the innovation process, 
### assuming that the scaled and mean adjusted residuals follow an AR(1) process.
### The function receives residuals for any given compartment/state and returns the resampled trajectory for that compartment/state.
resi_adj_boot = function(res, est_del_st, samp_ind){
  resi.adj = res - mean(res)
  s = sqrt(pmax(0,lowess((resi.adj)^2,f = 1/16)$y))
  # n = length(resi.adj)
  # ref = c(resi.adj[n:1], resi.adj, resi.adj[n:1])
  # s2 = sqrt(pmax(0,lowess(ref^2, f = bw_lowess/3)$y)[(n+1):(2*n)])
  # scale_factor = median(s2)/2
  # resi.adj.scale = resi.adj/(s2+scale_factor)

  resi.adj.scale = resi.adj/(s+.01)
  zeta = 0
  n = length(resi.adj.scale)
  for(t in 2:n){
    zeta_numer =+ (resi.adj.scale[t]*resi.adj.scale[t-1])
    zeta_denom =+ resi.adj.scale[t]^2
  }
  zeta = zeta_numer/ zeta_denom
  zeta = ifelse( zeta <= .98, zeta, .98)
  zeta = ifelse( zeta >= -.98, zeta, -.98)
  hat.resi.adj.scale = c(resi.adj.scale[1],rep(0,(n-1)))
  for(t in 2:n){
    hat.resi.adj.scale[t] = zeta*resi.adj.scale[t-1]
  }
  resi_innov = resi.adj.scale - hat.resi.adj.scale
  resi_innov.samp = resi_innov[samp_ind]
  resi.adj.scale.boot = resi_innov.samp + hat.resi.adj.scale
  #resi.adj_boot = (s2+scale_factor)*resi.adj.scale.boot
  resi.adj_boot = (s+.01)*resi.adj.scale.boot
  boot_samp = mean(res) + resi.adj_boot + est_del_st
  return(boot_samp)
}



### This function computes the bootstrapped version of all the compartments by calling the function resi_adj_boot  above.
resi_bootstrap = function(observed, estimated, Time){

  Max.pop = unique(observed$pop19)

  resi.Ct = diff(sqrt(abs(observed$C)))[1:(Time-2)] - diff(sqrt(abs(estimated$C)))
  resi.Rt = diff(sqrt(abs(observed$R)))[1:(Time-2)] - diff(sqrt(abs(estimated$R.Reported)))
  resi.Dt = diff(sqrt(abs(observed$D)))[1:(Time-2)] - diff(sqrt(abs(estimated$D)))
  resi.Ht = diff(sqrt(abs(observed$H)))[1:(Time-2)] - diff(sqrt(abs(estimated$H))) ##A
  resi.Qt = diff(sqrt(abs(observed$Q)))[1:(Time-2)] - diff(sqrt(abs(estimated$Q))) ##A

  samp_ind = sample(x = 1:(length(resi.Ct)), size = (length(resi.Ct)), replace = TRUE)
 
  delCt.adj.boot = resi_adj_boot(resi.Ct, diff(sqrt(abs(estimated$C))), samp_ind)
  delDt.adj.boot = resi_adj_boot(resi.Dt, diff(sqrt(abs(estimated$D))), samp_ind)
  delRt.adj.boot = resi_adj_boot(resi.Rt, diff(sqrt(abs(estimated$R.Reported))), samp_ind)
  delHt.adj.boot = resi_adj_boot(resi.Ht, diff(sqrt(abs(estimated$H))), samp_ind)
  delQt.adj.boot = resi_adj_boot(resi.Qt, diff(sqrt(abs(estimated$Q))), samp_ind)


  A = abs(estimated$A)
  C = c(abs(observed$C[1]),(cumsum(delCt.adj.boot) + sqrt(abs(observed$C[1])))^2)
  R = c(abs(observed$R[1]),(cumsum(delRt.adj.boot) + sqrt(abs(observed$R[1])))^2)
  D = c(abs(observed$D[1]),(cumsum(delDt.adj.boot) + sqrt(abs(observed$D[1])))^2)
  H = c(abs(observed$H[1]),(cumsum(delHt.adj.boot) + sqrt(abs(observed$H[1])))^2)
  Q = c(abs(observed$Q[1]),(cumsum(delQt.adj.boot) + sqrt(abs(observed$Q[1])))^2)
  RA = abs(estimated$R)[1:(Time-1)] - abs(estimated$R.Reported)[1:(Time-1)]
  S = pmax((rep(Max.pop,(Time-1)) - A[1:(Time-1)] - R - RA -Q[1:(Time-1)] -H[1:(Time-1)] -D -C),1)




  delC = diff(C)
  delR = diff(R)
  delH = diff(H)
  delD = diff(D)
  delQ = diff(Q)
  delS = diff(S)

  delTests = lowess((diff(observed$Test)), f = bw_lowess)$y
  TtoH_deno = lowess(H, f = bw_lowess)$y
  n_eff = min(length(delTests), length(TtoH_deno))
  TtoH = delTests[1:n_eff]/TtoH_deno[1:n_eff]

  var = list(delC = delC, delR = delR, delH = delH, delD = delD,
             A = A,
             Q = Q, H = H, TtoH = TtoH,
             S = S,
             R = R, D = D,
             C = C, Test = observed$Test[1:(Time -1)],
             pop19 = observed$pop19[1:(Time -1)],
             kappa= observed$kappa[1:(Time -1)])


  return(list(var = var, theta = estimated$theta.t, n = (Time -2)))
}


##For a given state of the US, we compute the boostrapped version of the trajectories of different compartments (given as a output of resi_bootstrap),
### and estimate the relevant parameters and epi-markers (by calling estimation_module function from the sourced script "est_module.R")
ncore <- 25
Nsim = 1000
###read all the variables
st = "Utah"
data.st = data_st(st) ##observed data
observed_states = data.st
load(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
estimated = est_st$est_states ## estimated data from estimation_module function in the sourced script "est_module.R"

#true parameters for the given state of the US
true.para.st <- list(gamma.true = est_st$est_para$est$gamma.hat,
                     rhoA.true = est_st$est_para$est$rhoA.hat,
                     rhoH.t.true = est_st$est_para$est$rhoH.final,
                     alpha.true = (mean(est_st$est_para$est$sqrt.alpha.kappa.t/data.st$var$kappa[1:length(est_st$est_para$est$sqrt.alpha.kappa.t)]))^2,
                     delta.t.true = est_st$est_para$mortality$mort_lowess_rate,
                     phi.t.true = est_st$est_para$est$phi.hat.smooth)


#The output containing the bootstrapped "observed" compartments and the corresponding estimated parameters/ compartments
est.list <- mclapply(1:Nsim, function(i){
  boot_sample = resi_bootstrap( observed = observed_states$var,
                                estimated = estimated,
                                Time = data.st$n)


  boot_est = estimation_module(dat = boot_sample,
                               numBlock = est_st$tuning$numBlock ,
                               wd.length = est_st$tuning$wd.length,
                               bw_lowess = bw_lowess,
                               boot = TRUE)

  return(list(boot_sample = boot_sample, boot_est = boot_est))
}, mc.cores = ncore)

File <- sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_bootstrap_A.Rda", st, st, st)
if (file.exists(File)) stop(File, " already exists")
dir.create(dirname(File), showWarnings = FALSE)
save(est.list, file = 
       sprintf("../Results for USA States/%s/%s_updated_boot_final/%s_bootstrap_A.Rda", st, st, st))

