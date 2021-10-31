




rm(list=ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(parallel)
load("./set1_wd_30.RData")
source("./func_101221.R")
source('./fitting_function_101620.R')

## dat.sim.new is obtained from simulate_process.R
## dat.sim.new.est is obtained from fit_data.R

Nsim <- 1000
ncore <- 15

est.BS.list <- mclapply(1:Nsim, function(i){
  
  set.seed(456*i)
  
  dat.sim.new.bs <- lapply(1:10, function(i){
    return(resi_bootstrap(
      observed = dat.sim.new[[i]]$data,
      estimated =  dat.sim.new.est$est_states[[i]]$data, 
      theta = dat.sim.new.est$est_states[[i]]$theta,
      Time.use = 94, Max.pop = Max.pop[i]))
  })
  
  
  boot_est = estimation_module_agg(dat = dat.sim.new.bs,
                                   numBlock = 9,
                                   wd.length = 30,Time=100) 
  #print(true.para)
  return(list(boot_sample = dat.sim.new.bs, boot_est = boot_est))
}, mc.cores = ncore)

results <- list(est.BS.list = est.BS.list,
                dat.sim.new = dat.sim.new,
                dat.sim.new.est = dat.sim.new.est)


save(results, file='agg_res_bs_101321_poisson.Rda')


#save(results, file='agg_res_bs_101321_gumbel.Rda')