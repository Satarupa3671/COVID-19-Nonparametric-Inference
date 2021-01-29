
###Based on the output of the resi_boot.R, computes the mean estimates for different parameters across different states. 
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
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


#source("./est_module_MN.R")
source("./est_module.R")

load("../USA_data_processing/covid_df.Rda")
st_list = c("Arizona", "Arkansas",  "Delaware", "Idaho", "Iowa", #Minnesota,
            "Nebraska","Ohio", "Oklahoma", "Pennsylvania",
             "South Dakota", "Tennessee", "Texas", "Utah", "Wisconsin")
estd_data = list()
for(st in st_list){ 
  #c("Utah", "Iowa", "Massachusetts","Pennsylvania","Idaho","North Carolina")){ #
   # "Minnesota"){
  load(sprintf("../Results for USA States/%s/%s_est_st.Rda", st, st))
  estd_data[[st]] = est_st
}

gamma_est = sapply(estd_data, function(st){
  st$est_para$est$gamma.hat
})

rhoA_est = sapply(estd_data, function(st){
  st$est_para$est$rhoA.hat
})


alpha_est = sapply(names(estd_data), function(st_name){
  n = length(estd_data[[st_name]]$est_para$est$sqrt.alpha.kappa.t)
  round(mean(((estd_data[[st_name]]$est_para$est$sqrt.alpha.kappa.t) / (data_st(st_name)$var$kappa[1:n]))^2),4)
})


delta_est_mean =  sapply(estd_data, function(st){
  round(mean(st$est_para$mortality$mort_lowess_rate),4)
})


rhoH_est_mean = sapply(estd_data, function(st){
  round(mean(st$est_para$est$rhoH.hat.smooth),4)
})

phi_est_mean = sapply(estd_data, function(st){
  round(mean(st$est_para$est$phi.hat.smooth),4)
})


theta_est_mean =  sapply(estd_data, function(st){
  round(mean(st$est_para$est$theta.hat.t),4)
})

table_df = as.data.frame(cbind(st_list, gamma_est, rhoA_est , alpha_est, delta_est_mean, rhoH_est_mean , phi_est_mean, theta_est_mean))

print(xtable(table_df))

