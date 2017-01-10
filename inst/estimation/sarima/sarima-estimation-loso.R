## regional SARIMA fits
## started: 5 Oct 2016
## Nicholas Reich

library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(awes)

## 5 Oct 2016: 
## need to install forecast from GitHub to take advantage of recent bug fix.
## devtools::install_github("robjhyndman/forecast")
library(forecast)

### Load data set and set variables describing how the fit is performed
## Load data for nationally reported influenza like illness
data <- read.csv("data-raw/allflu-cleaned.csv", stringsAsFactors = FALSE)
data$time <- as.POSIXct(data$time)

## where to store SARIMA fits
path <- paste0("inst/estimation/sarima/fits/")

## do all regional fits in parallel
library(doMC)
registerDoMC(4)
foreach(reg_num=c("X", 1:10)) %dopar%
  fit_region_sarima(data,
    reg_num = reg_num,
    first_test_season = "2011/2012", # 5 test seasons: 2011/2012, 2012/2013, 2013/2014, 2014/2015, 2015/2016
    seasonal_difference = TRUE,
    transformation = "log",
    path = path)
