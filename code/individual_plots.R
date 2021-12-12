### runs the estimation module for a few states and saves the respective plots. 

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(scales)
library(lemon)

source("./fitting_functions.R")
source("./est_plot_func.R")
source("./est_module.R")
load("../USA_data_processing/covid_df.Rda")
bw_lowess = 1/16

##data_st is the function to read the data for a particular state "st" and return the smoothed trajectories of (C,D,H,Q,R,DelC,DelH,DelR,DelQ,DelD,Test,Kappa) as a list

##opt_wd contains the optimal wd.length and numBlock for a select few states.
opt_wd = read.csv("../USA_data_processing/opt_wd_list.csv",header = TRUE)
#opt_wd = read.csv("../USA_data_processing/opt_wd_list2.csv",header = TRUE)

for (st in opt_wd$states){
  numBlock = opt_wd$numBl_st[opt_wd$states == st]
  wd.length = opt_wd$wdL_st[opt_wd$states == st]
  plot_for_estimation(st, numBlock, wd.length)
}
######
source("./est_module_MN.R")
st = "Minnesota"
bw_lowess = 1/16
plot_for_estimation(st, numBlock, wd.length)
