rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magrittr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(lubridate)
library(plyr)
library(gridExtra)

###Total Cases for each state.

df1 = read.csv("~/Dropbox/covid/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv",header=FALSE)
dates=unlist(df1[1,setdiff(1:(dim(df1)[2]),1:11)], use.names=FALSE)
df2 = read.csv("~/Dropbox/covid/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv",header=FALSE, skip = 1)
df_dates = df2[-c(1:6,8:11)]
colnames(df_dates)[1] = "state"
colnames(df_dates)[2: dim(df_dates)[2]] = as.character(dates)
unique(df_dates$state)

data_spl = lapply(unique(df_dates$state), function(prov){
  subset(df_dates, df_dates$state == prov)
})

cases = t(sapply(1: length(data_spl), function(prov){
  as.numeric(colSums(data_spl[[prov]][,-1]))
}))


cases_df = as.data.frame(cbind(as.character(unique(df_dates$state)),cases))
colnames(cases_df) = c("state",as.character(dates))

#put it into long format
cases_df1 = cases_df %>% melt(id.vars = "state")
cases_df1 = cases_df1[order(cases_df1$state),]
names(cases_df1) = c("state", "date", "total_cases")
cases_df1$date = cases_df1$date %>% as.Date("%m/%d/%y")
cases_df1$total_cases = cases_df1$total_cases%>% as.numeric()
cases_df1 = cases_df1 %>% mutate(new_cases = total_cases - dplyr::lag(total_cases)) #Create new deaths


###Total Deaths for each state. 

df1 = read.csv("~/Dropbox/covid/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv",header=FALSE)
dates=unlist(df1[1,setdiff(1:(dim(df1)[2]),1:11)], use.names=FALSE)
df2 = read.csv("~/Dropbox/covid/data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv",header=FALSE, skip = 1)
df_dates = df2[-c(1:6,8:11)]
colnames(df_dates)[1] = "state"
colnames(df_dates)[2: dim(df_dates)[2]] = as.character(dates)


data_spl = lapply(unique(df_dates$state), function(prov){
  subset(df_dates, df_dates$state == prov)
})

deaths = t(sapply(1: length(data_spl), function(prov){
  as.numeric(colSums(data_spl[[prov]][,-1]))
}))


deaths_df = as.data.frame(cbind(as.character(unique(df_dates$state)),deaths))
colnames(deaths_df) = c("state",as.character(dates))

#put it into long format
deaths_df1 = deaths_df %>% melt(id.vars = "state")
deaths_df1 = deaths_df1[order(deaths_df1$state),]
names(deaths_df1) = c("state", "date", "total_deaths")
deaths_df1$date = deaths_df1$date %>% as.Date("%m/%d/%y")
deaths_df1$total_deaths = deaths_df1$total_deaths%>% as.numeric()
deaths_df1 = deaths_df1 %>% mutate(new_deaths = total_deaths - dplyr::lag(total_deaths)) #Create new deaths

##this dataframe consists of the total cases and deaths for each state.
covid_df = left_join(cases_df1, deaths_df1, by = c("state", "date"))
covid_df$state=as.vector(covid_df$state)


#Obtain state abbreviation to pull COVID tracking dataset
df_codes=read.csv("~/Dropbox/covid/data/EconomicTracker/data/GeoIDs - State.csv", header=TRUE)
df_codes=df_codes[,c("statename","stateabbrev","statefips","state_pop2019")]
colnames(df_codes)[which(colnames(df_codes)=="statename")]="state"
colnames(df_codes)[which(colnames(df_codes)=="state_pop2019")]="pop19"
colnames(df_codes)[which(colnames(df_codes)=="stateabbrev")]="stateabbrv"
df_codes$state=as.character(df_codes$state)
df_codes$stateabbrv=as.character(df_codes$stateabbrv)
df_codes$state[df_codes$state == "District Of Columbia"] = "District of Columbia"
##This dataframe contains additional information on state abbv.
covid_df1 = inner_join(covid_df, df_codes, by = c("state"))


#Read tracking project data
setwd("~/Dropbox/covid/data/covid_tracking/states_data/")
filelist = list.files("~/Dropbox/covid/data/covid_tracking/states_data/", pattern = "*.csv")
tracking_state = do.call(rbind,
                         lapply(filelist, read.csv, header = TRUE))


tracking_state = tracking_state %>% 
  select(c("date","state","positive","negative","hospitalizedCurrently",
           "hospitalizedCumulative","recovered"))

tracking_state = tracking_state %>% 
  mutate(date = as.Date(as.character(tracking_state$date)))

colnames(tracking_state)=c("date","stateabbrv","positive","negative","hospitalizedCurrently",
                           "hospitalizedCumulative","recovered")
tracking_state$stateabbrv=as.character(tracking_state$stateabbrv)
tracking_state$total_tests = tracking_state$positive + tracking_state$negative

##This dataframe contains additional information on hosp.,rec.,tests.

covid_df2 = covid_df1 %>% inner_join(tracking_state, by = c("stateabbrv","date"))

#Read Google mobility data
df_google = read.csv("~/Dropbox/covid/data/EconomicTracker/data/Google Mobility - State - Daily.csv",header = TRUE)
df_google=df_google %>% na_if(".")
df_google[,5:ncol(df_google)]=apply(df_google[,5:ncol(df_google)], 2, as.numeric)
df_google$date=paste(paste(df_google$year, df_google$month, sep = "-"), df_google$day, sep = "-") %>% as.Date("%Y-%m-%d")
df_google=df_google[,-c(1,2,3)]

df_google$average <- rowMeans(df_google[,2:7], na.rm = T)

## This dataframe has further information on google mobility data
covid_final = inner_join(covid_df2,df_google, by = c("date","statefips"))
covid_final = covid_final[order(covid_final$state),]

#View(covid_df[covid_df$state == "New York",])

save(covid_final,file="~/Dropbox/covid/data/USA_data_processing/covid_df.Rda")


