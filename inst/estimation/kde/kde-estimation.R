## KDE estimation
## 10/7/2016 - Nicholas Reich - file created
## 1/10/2017 - Evan Ray - moved to awes repository, minor modifications
## 1/11/2017 - Nicholas Reich - modified to create loso fits/predictions

library(awes)

fludat <- read.csv('data-raw/allflu-cleaned.csv')

region_strings <- c("X", paste("Region", 1:10))

## fit LOSO models on training seasona
library(doMC)
registerDoMC(4)
foreach(reg=region_strings) %dopar%
    fit_region_kdes(fludat, region=reg,
                    first_test_year = 2011, first_test_week = 31,
                    path = "inst/estimation/kde/fits/")

## make leave-one-season-out predictions for training seasons
n_sims <- 10000
reg=region_strings[1]
#foreach(reg=region_strings) %dopar%
      for(reg in region_strings){
    predict_region_kde(fludat, region=reg, 
                       path = "inst/estimation/kde/fits/",
                       n_sim = n_sims)
 }

source("inst/estimation/kde/check-kde-predictions.R")