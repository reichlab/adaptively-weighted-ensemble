## KDE estimation
## 10/7/2016 - Nicholas Reich - file created
## 1/10/2017 - Evan Ray - moved to awes repository, minor modifications
## 1/11/2017 - Nicholas Reich - modified to create loso fits/predictions
## 1/16/2017 - Nicholas Reich - modified for standardized LOSO prediction function
## 1/18/2017 - Nicholas Reich - random seeds and re-organizing of code

## set to root directory of package code
setwd("~/Documents/code_versioned/adaptively-weighted-ensemble/")

library(awes)
library(rstream)

## setup
fludat <- read.csv('data-raw/allflu-cleaned.csv')
region_strings <- c("X", paste("Region", 1:10))
fit_path <- "inst/estimation/kde/fits/"
method <- "kde"
analysis_time_seasons <- paste0(1997:2010, "/", 1998:2011)

## setup for random seeds
rngstream <- get_initial_rng_substream()

library(doMC)
registerDoMC(11)
foreach(reg=region_strings) %dopar% {
    
    ## fit LOSO models on training seasona
    # reg = region_strings[1]
    fit_region_kdes(fludat, region=reg,
                    first_test_year = 2011, first_test_week = 31,
                    path = fit_path)
    
    ## make leave-one-season-out predictions for training seasons
    # ## for debugging
    # debug(get_log_scores_via_direct_simulation)
    # region = region_strings[1]
    # analysis_time_season = "2000/2001"
    for(analysis_time_season in analysis_time_seasons) {
        ## set random seed from substream
        get_rng_substream(
            rngstream = rngstream,
            method = method,
            region = reg,
            season = analysis_time_season,
            set_rng = TRUE)
        
        ## get and save log scores
        get_log_scores_via_direct_simulation(
            analysis_time_season = analysis_time_season,
            first_analysis_time_season_week = 10, # == week 40 of year
            last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
            region = reg,
            prediction_target_var = "weighted_ili",
            incidence_bins = data.frame(
                lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
                upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
            incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
            n_sims = 100000,
            model_name = "kde",
            fits_path = fit_path,
            prediction_save_path = "inst/estimation/loso-predictions/"
        )
    }
}
    
source("inst/estimation/kde/check-kde-predictions.R")
