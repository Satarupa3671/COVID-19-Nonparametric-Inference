rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)


Global_Mobility_Report = read.csv("../../data/USA_data_processing/Global_Mobility_Report.csv")

mobility <- Global_Mobility_Report %>%
  filter(country_region_code=='US') %>%
  filter(!is.na(iso_3166_2_code)) %>%
  select(-sub_region_2,-metro_area,-census_fips_code) 


mobility$average <- rowMeans(mobility[,6:11], na.rm = T)


plot(1:164,mobility$average[1:164])

write.csv(mobility,file = "~/Dropbox/covid/data/USA_data_processing/US_mobility_021520_081720.csv")

