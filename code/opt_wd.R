###calculates the optimal wd.length and numBlock parameters from the observed data in view of the most updated estimation code. (We only use the optimal wd.length and numBlock obtained from the results of this file, and don't need it otherwise)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
#library(tidyverse)
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
bw_lowess = 1/16

##Finds the optimal numBlock and wd.length parameters for fitting the model without using the kappa (mobility) information.
opt_numblock = function(dat,numBlock, wd.length){

  est_st = estimation_module(dat = dat, numBlock = numBlock, wd.length = wd.length , boot = FALSE)
  
  est1 = est_st$est_para$est
  s1 =  (mean((lowess(diff(est1$A.hat.t), f = bw_lowess)$y - lowess(est1$delA.hat.t,f = bw_lowess)$y)^2))/
    (abs(mean(lowess(diff(est1$A.hat.t), f = bw_lowess)$y) + mean(lowess(est1$delA.hat.t,f = bw_lowess)$y)))
  
  s2 =  (mean((lowess(est1$new.infect.t,f = bw_lowess)$y - lowess(dat$var$delC[2:est1$n.t],f = bw_lowess)$y)^2))/
    (abs(mean(lowess(est1$new.infect.t,f = bw_lowess)$y) + mean(lowess(dat$var$delC[2:est1$n.t],f = bw_lowess)$y)))
  
  s3 =  (mean((lowess(cumsum(est1$new.infect.t), f = bw_lowess)$y - lowess(cumsum(dat$var$delC[2:est1$n.t]),f = bw_lowess)$y)^2))/
    (abs(mean(lowess(cumsum(est1$new.infect.t), f = bw_lowess)$y) + mean(lowess(cumsum(dat$var$delC[2:est1$n.t]),f = bw_lowess)$y)))
  
  return(s1 + s2 + s3)
}

st_list = c("Arizona", "Arkansas",  "Delaware", "Idaho", "Iowa", #Minnesota,
            "Nebraska","Ohio", "Oklahoma", "Pennsylvania",
            "South Dakota", "Tennessee", "Texas", "Utah", "Wisconsin")


for (st in st_list){
data.st = data_st(st)
numBlock.grid = 10:14
wd.length.grid = seq(15,60,1)
ncore = 5
eval.wd.length = sapply(numBlock.grid , function(numbl){
  mclapply(wd.length.grid, function(wd){
    opt_numblock(data.st,numbl,wd)
  },mc.cores = ncore)
})

eval.wd.length2 = matrix(c(unlist(eval.wd.length)), 
                         nrow= length(wd.length.grid), ncol = length(numBlock.grid),byrow = F)
pos = which(eval.wd.length2 == min(eval.wd.length2,na.rm = TRUE), arr.ind=TRUE)
print(paste(st, "--> ", wd.length.grid[pos[1]],numBlock.grid[pos[2]]))  
}



