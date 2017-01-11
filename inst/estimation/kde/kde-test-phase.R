## KDE testing phase models and predictions
## 1/11/2017 - Nicholas Reich - created, adapted from KoTflu20162017 package repository

## create set of testing phase fits to use in generating predictions, see issue #17
library(doMC)
registerDoMC(4)
foreach(reg=region_strings) %dopar%
    fit_region_kdes(fludat, region=reg,
                    first_test_year = 2016, first_test_week = 31, 
                    ## setting first_test_week to 31 ensures that we create a fit for 2016-2017
                    path = "inst/estimation/kde/loso-fits/")

## move fits for current year to current-fits folder
mv_command <- paste0("mv inst/estimation/kde/loso-fits/*2016-2017.rds inst/estimation/kde/current-fits/")
system(mv_command)

