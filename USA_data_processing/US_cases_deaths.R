rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)

###Total Cases for each state. We eliminate 6 non-populous states/cruise ships.

df1 = read.csv("../../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv",header=FALSE)
dates=unlist(df1[1,setdiff(1:(dim(df1)[2]),1:11)], use.names=FALSE)
df2 = read.csv("../../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv",header=FALSE, skip = 1)
df_dates = df2[-c(1:6,8:11)]
colnames(df_dates)[1] = "province"
colnames(df_dates)[2: dim(df_dates)[2]] = as.character(dates)
unique(df_dates$province)

drops = c("American Samoa","Guam","Northern Mariana Islands","Virgin Islands","Diamond Princess","Grand Princess")
df_dates = subset(df_dates, !(province %in% drops))
data_spl = lapply(unique(df_dates$province), function(prov){
  subset(df_dates, df_dates$province == prov)
})

cases = t(sapply(1: length(data_spl), function(prov){
  as.numeric(colSums(data_spl[[prov]][,-1]))
}))


cases_df = as.data.frame(cbind(as.character(unique(df_dates$province)),cases))
colnames(cases_df) = c("province",as.character(dates))

#put it into long format
cases_df1 = cases_df %>% melt(id.vars = "province")
cases_df1 = cases_df1[order(cases_df1$province),]
names(cases_df1) = c("province", "date", "total_cases")
cases_df1$date = cases_df1$date %>% as.Date("%m/%d/%y")



###Total Deaths for each state. We eliminate 6 non-populous states/cruise ships.

df1 = read.csv("../../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv",header=FALSE)
dates=unlist(df1[1,setdiff(1:(dim(df1)[2]),1:11)], use.names=FALSE)
df2 = read.csv("../../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv",header=FALSE, skip = 1)
df_dates = df2[-c(1:6,8:11)]
colnames(df_dates)[1] = "province"
colnames(df_dates)[2: dim(df_dates)[2]] = as.character(dates)


drops = c("American Samoa","Guam","Northern Mariana Islands","Virgin Islands","Diamond Princess","Grand Princess")
df_dates = subset(df_dates, !(province %in% drops))
data_spl = lapply(unique(df_dates$province), function(prov){
  subset(df_dates, df_dates$province == prov)
})

deaths = t(sapply(1: length(data_spl), function(prov){
  as.numeric(colSums(data_spl[[prov]][,-1]))
}))


deaths_df = as.data.frame(cbind(as.character(unique(df_dates$province)),deaths))
colnames(deaths_df) = c("province",as.character(dates))

#put it into long format
deaths_df1 = deaths_df %>% melt(id.vars = "province")
deaths_df1 = deaths_df1[order(deaths_df1$province),]
names(deaths_df1) = c("province", "date", "total_deaths")
deaths_df1$date = deaths_df1$date %>% as.Date("%m/%d/%y")

covid_df = left_join(cases_df1, deaths_df1, by = c("province", "date"))
covid_df$province=as.vector(covid_df$province)
covid_df$total_cases = covid_df$total_cases %>% as.numeric()
covid_df$total_deaths = covid_df$total_deaths %>% as.numeric()


newcases = function(df){
  df %>%
    mutate(new_cases = total_cases - dplyr::lag(total_cases),
           new_deaths = total_deaths - dplyr::lag(total_deaths))
}

start_in_april = function(df){ 
  df %>% filter(date >= "2020-04-13")
}


make_time = function(df){
  if(nrow(df)>0){
    maxT = max(nrow(df),1)
    df %>% mutate(t = (1:maxT)-1)
  }
}

process = function(df){ #Receives a country for example: \code{covid_df %>% filter(country == "Chile")}
  df %>% newcases %>% 
    start_in_april %>% 
    make_time
}

covid_df = ddply(.data = covid_df, .(province), process) 

write.csv(covid_df,file = "~/Dropbox/covid/data/USA_data_processing//US_cases_deaths.csv")
