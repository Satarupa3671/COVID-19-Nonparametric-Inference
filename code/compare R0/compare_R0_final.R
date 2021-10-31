### Function to generate new R0 plots corresponding to the updated data ("covid_df_2020March7_to_2021March7.Rda")
### This plot is supplied in the reference letter to answer a reviewer.
### Useful functions are 


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)
library(scales)
library(lemon)
library(plyr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(gridExtra)


source("../fitting_functions.R")
source("./est_plot_func_R0.R")
source("./est_module_R0.R")

load("../../USA_data_processing/covid_df.Rda")
covid_old = covid_final
load("../../USA_data_processing/covid_df_2020March7_to_2021March7.Rda")
covid_new = covid_final

load("../../Results for USA States/Utah/Utah_est_st.Rda")
est_st_old = est_st

st = "Utah"
covid.state = subset(covid_new, state == st)
covid.state_old = subset(covid_old, state == st)

bw_lowess_new = .033
opt_wd = read.csv("../../USA_data_processing/opt_wd_list1.csv",header = TRUE)

calc_R0 = function(st, bw_lowess){
  data.st = data_st(st,bw_lowess)
  numBlock = opt_wd$numBl_st[opt_wd$states == st]
  wd.length = opt_wd$wdL_st[opt_wd$states == st]
  est_st = estimation_module(dat = data.st, numBlock = numBlock , wd.length = wd.length, bw_lowess = bw_lowess, boot = FALSE)
  print(paste0("gamma = ", est_st$est_para$est$gamma.hat))
  print(paste0("rhoA = ", est_st$est_para$est$rhoA.hat))
  print("\n")
  R0_numer = est_st$est_para$est$sqrt.alpha.kappa.t^2 * unique(data.st$var$pop19)
  
  R0_denom = (est_st$est_states$S + est_st$est_states$A + est_st$est_states$R) * 
    (est_st$est_para$est$theta.hat.t + est_st$est_para$est$gamma.hat + est_st$est_para$est$rhoA.hat)
  R0 = R0_numer/R0_denom[1:length(R0_numer)]
  return(R0)
}
data.st = data_st("Utah", 1/16)

##R0 old
R0_numer_old = est_st_old$est_para$est$sqrt.alpha.kappa.t^2 * unique(covid.state_old$pop19)

R0_denom_old = (est_st_old$est_states$S + est_st_old$est_states$A + est_st_old$est_states$R) * 
  (est_st_old$est_para$est$theta.hat.t + est_st_old$est_para$est$gamma.hat + est_st_old$est_para$est$rhoA.hat)
R0_old = R0_numer_old/R0_denom_old[1:length(R0_numer_old)]

##R0_new
R0_new_bw = calc_R0("Utah", bw_lowess_new)
##combining it all
R0_combined = data.frame (Days = 1:length(R0_new_bw), R0_new = R0_new_bw)
R0_combined$dates1 = as.Date(data.st$dates[1:length(R0_new_bw)], format = "%m/%d/%Y")

n_extra1 = length(R0_new_bw) - length(R0_old)
R0_combined$R0_old = c(R0_old, rep(NA,n_extra1))

##ihme

#ihme = read.csv("../data/ihme.csv",header = TRUE)
# ihme_Utah = (ihme[ihme$location_name == "Utah",]$mean)
# ihme_Utah = ihme_Utah[!is.na(ihme_Utah)]
ihme = read.csv("../../data/ihme.csv",header = TRUE)
ihme_Utah = ihme[ihme$location_name == "Utah",]
ihme_Utah$date = as.Date(as.character(ihme_Utah$date))#, "%m/%d/%Y")

ihme_Utah$date = ihme_Utah$date[1:which.max(ihme_Utah$date) ]

c(min(R0_combined$dates1), max(R0_combined$dates1)) ##"2020-05-07" ---"2021-03-05"
c(min(ihme_Utah$date),max(ihme_Utah$date)) ##"2019-12-12"-- "2021-02-28"
dates_effective = c(max(min(R0_combined$dates1), min(ihme_Utah$date)),
                    max(max(R0_combined$dates1), max(ihme_Utah$date))) ##"2020-05-07" ----"2021-03-05"
ihme_Utah = ihme_Utah[ihme_Utah$date >= dates_effective[1] & ihme_Utah$date <= dates_effective[2],]
R0_ihme = ihme_Utah$mean

R0_combined2 = R0_combined[R0_combined$dates1 >= dates_effective[1] & R0_combined$dates1 <= dates_effective[2],]

R0_combined2$R0_ihme =  c(R0_ihme, rep(NA, 5))
###
rt_live = read.csv("../../data/rt.csv", header = TRUE)
Utah_R0 = as.data.frame(rt_live[rt_live$region == "UT",])


R0_df2 = data.frame (Days = 1:length(Utah_R0$median), Rt_live = Utah_R0$median)
R0_df2$dates1 = as.Date(Utah_R0$date)[1:length(Utah_R0$median)]
c(min(R0_df2$dates1), max(R0_df2$dates1)) ##dates they have data on "2020-02-29" ---"2021-01-26"
c(min(R0_combined2$dates1), max(R0_combined2$dates1))##dates we have data on "2020-05-07"--- "2021-03-05"

R0_combined2$dates1[1] < R0_df2$dates1[1] #FALSE
R0_combined2$dates1[length(R0_combined2$dates1)] > R0_df2$dates1[length(R0_df2$dates1)] ##TRUE

#They have data from Feb 29 to Jan 26. We have only from May 7  to Dec 4. We need to subset only that time.

dates_effective = c(min(R0_combined2$dates1),max(R0_combined2$dates1)) ##"2020-05-07" ----"2021-03-05"
R0_df2 = R0_df2[R0_df2$dates1>= dates_effective[1] & R0_df2$dates1<= dates_effective[2],]

R0_df2 = R0_df2[(R0_df2$dates1 >= R0_combined2$dates1[1]) & (R0_df2$dates1 <= R0_combined2$dates1[length( R0_combined2$dates1)]),]
extra_dates1 = seq(max(R0_df2$dates1)+1, dates_effective[2],1)
extra_Days = sapply(1:length(extra_dates1), function(dt) max(R0_df2$Days)+dt)
extra_Rt_live = rep(NA,length(extra_Days))
extra_df_R0_df2 = data.frame(Days = extra_Days, Rt_live = extra_Rt_live, dates1 = extra_dates1)
R0_combined3 = rbind(R0_df2,extra_df_R0_df2)
R0_combined3$R0_corresponding =  R0_combined2$R0_new
R0_combined3$R0_old =  R0_combined2$R0_old
R0_combined3$R0_ihme =  R0_combined2$R0_ihme






setwd("../../data/r_values/")
filelist = list.files(,pattern = "*r_values_us.csv")

filelist1 =  lapply(filelist, read.csv, header = TRUE)
pos = min( which(sapply(1: length(filelist1), function(ind) (ncol(filelist1[[ind]]) > 4 )) == TRUE))
filelist2 =do.call(rbind, lapply(pos :length(filelist1), function(ind){
  cbind(as.character(filelist1[[ind]]$region),
        as.numeric(as.character(filelist1[[ind]]$initial_r_0)),
        as.numeric(as.character(filelist1[[ind]]$post_mitigation_r)),
        as.character(filelist1[[ind]]$cur_date))
}))
filelist2 = as.data.frame(filelist2)
colnames(filelist2) = c('region','initial_r_0', 'current_r', 'cur_date')

Utah_R0 = as.data.frame(filelist2[filelist2[,1] == "UT",])
Utah_R0$initial_r_0 = as.numeric(as.character(Utah_R0$initial_r_0))
Utah_R0$current_r = as.numeric(as.character(Utah_R0$current_r))
Utah_R0$cur_date = as.Date(Utah_R0$cur_date)
Utah_R0 = Utah_R0[-nrow(Utah_R0),]



c(min(R0_combined$dates1), max(R0_combined$dates1)) ##"2020-05-07" ---"2021-03-05"
c(min(Utah_R0$cur_date),max(Utah_R0$cur_date)) ##"2020-05-24"--- "2020-10-04"

extra_dates1 = seq(dates_effective[1],min(Utah_R0$cur_date)-1,1)
extra_dates2 = seq(max(Utah_R0$cur_date)+1, dates_effective[2],1)
c(extra_dates1,Utah_R0$cur_date, extra_dates2)
R0_combined3$R0_SIER =  c(rep(NA,length(extra_dates1)),Utah_R0$current_r,rep(NA,length(extra_dates2)))

break.vec <- c(R0_combined3$dates1[1],
               seq(from = R0_combined3$dates1[1] , to = R0_combined3$dates1[length(R0_combined3$dates1)],
                   by="month"), R0_combined3$dates1[length(R0_combined3$dates1)])

R0_plot_combined = ggplot(data = R0_combined3)+
  geom_line(mapping = aes(y = Rt_live, x = dates1), size = 1, col = 'blue', linetype = 2) +
  geom_line(mapping = aes(y = R0_SIER, x = dates1), size = 1, col = 'darkgreen', linetype = "3313") +
  geom_line(mapping = aes(y = R0_ihme, x = dates1), size = 1, col = 'magenta',linetype = 6) +
  geom_line(mapping = aes(y = R0_corresponding, x = dates1), size = 1, col = 'cyan', linetype = "solid") +
  #geom_line(mapping = aes(y = R0_old, x = dates1), size = 1, col = 'red', linetype = "dotted") +
  geom_hline(yintercept = 1, col ='black') +
  labs(x = "Time", y = expression(R[0])) +
  theme_bw() +
  scale_x_date(breaks = break.vec, #date_breaks("1 months"),
               labels = date_format("%b %d, %y")) +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
        axis.title = element_text(size = 26))+
  theme(legend.position = "none", #c(.5,.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 24))
R0_plot_combined
ggsave(file ="../../../Results for USA States/Utah/Utah_updated/R0_reviesd.eps")

####







break.vec <- c(R0_combined$dates1[1],
               seq(from = R0_combined$dates1[1] , to = R0_combined$dates1[length(R0_combined$dates1)], 
                   by="month"), R0_combined$dates1[length(R0_combined$dates1)])
R0_plot_combined = ggplot(data = R0_combined2)+
  geom_line(mapping = aes(y = R0_new, x = dates1), size = 1, col = 'red', linetype = "solid") +
  geom_line(mapping = aes(y = R0_old, x = dates1), size = 1, col = 'blue', linetype = "solid") +
  geom_line(mapping = aes(y = R0_ihme, x = dates1), size = 1, col = 'magenta', linetype = "solid") +
  geom_hline(yintercept = 1, col ='black') +
  labs(x = "Time", y = expression(R[0])) +
  theme_bw() +
  scale_x_date(breaks = break.vec, #date_breaks("1 months"),
               labels = date_format("%b %d, %y")) +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
        axis.title = element_text(size = 26))+
  theme(legend.position = "none", #c(.5,.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 24))
R0_plot_combined
ggsave(file ="../Results for USA States/Utah/Utah_updated/R0_overlaid2.pdf")

