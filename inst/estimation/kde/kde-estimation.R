## KDE estimation
## 10/7/2016 - Nicholas Reich - file created
## 1/10/2017 - Evan Ray - moved to awes repository, minor modifications

library(awes)

fludat <- read.csv('data-raw/allflu-cleaned.csv')

region_strings <- c("X", paste("Region", 1:10))

library(doMC)
registerDoMC(4)
foreach(reg=region_strings) %dopar%
    fit_region_kdes(fludat, region=reg,
                    first_test_year = 2011, first_test_week = 31,
                    path = "inst/estimation/kde/fits/")
