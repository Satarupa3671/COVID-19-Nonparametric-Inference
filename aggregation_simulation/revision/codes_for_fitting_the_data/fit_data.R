## aggregation simulation -- fit the data

rm(list=ls())
library(tidyverse)
library(parallel)
source("./func_101221.R")
source('./fitting_function_101620.R')

## dat.sim.new is obtained from simulate_process.R

dat.sim.new.est = estimation_module_agg(dat = dat.sim.new, 
                                        numBlock = 9, 
                                        wd.length = 30,Time=100) 
