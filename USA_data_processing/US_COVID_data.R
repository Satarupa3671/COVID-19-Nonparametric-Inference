rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
lst = list.files("../../data/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports_us", pattern = ".csv")
lst1 = substr(lst, 1, nchar(lst)-6)
lst1 = gsub("-", "/", lst1)
dates = as.Date(lst1, "%m/%d/%y")

setwd("../../data/COVID-19/csse_covid_19_data/csse_covid_19_daily_reports_us")
dat = sapply(lst, read.csv, simplify = FALSE)

drops = c("American Samoa","Guam","Northern Mariana Islands","Virgin Islands","Diamond Princess","Grand Princess","Recovered")
dat_each_day = lapply(1:length(dat), function(day){
  df_dates = subset(dat[[day]], !(Province_State %in% drops))
  df_dates = df_dates[,-c(2:5,15,16)]
  df_dates = cbind(date = rep(dates[day],each = nrow(df_dates)), df_dates)
  
})
for(i in 2:length(dat_each_day)){
 names(dat_each_day[[i]]) <- names(dat_each_day[[1]]) 
}
dat_all_days = do.call(rbind,dat_each_day)
dat_all_days = dat_all_days[order(dat_all_days$Province_State),]


write.csv(dat_all_days,file = "~/Dropbox/covid/data/USA_data_processing/US_COVID_data.csv")

