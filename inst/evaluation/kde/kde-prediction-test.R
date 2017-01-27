## create testing phase predictions for KDE model
## 1/19/2017 - Nick Reich : created file
library(awes)
library(rstream)

## set to root directory of package code
setwd("~/Documents/code_versioned/adaptively-weighted-ensemble/")

## setup
region_strings <- c("X", paste("Region", 1:10))
fit_path <- "inst/estimation/kde/fits/"
prediction_path <- "inst/evaluation/test-predictions/"
method <- "kde"
analysis_time_season <- "2011/2012"

seasons_to_copy <- paste0(2012:2015, "-", 2013:2016)

## setup for random seeds
rngstream <- get_initial_rng_substream(seed=95381) ## seed generated on random.org

library(doMC)
registerDoMC(11)
foreach(reg=region_strings) %dopar% {
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
        prediction_save_path = prediction_path
    )
    
    ## since KDE predictions are only based on historical data, we can copy predictions from the first season
    for(season in seasons_to_copy){
        if(identical(reg, "X")) {
          reg <- "National"
        }
        orig_file <- paste0(prediction_path, method, "-", gsub(" ", "", reg), "-2011-2012-loso-predictions.rds")
        new_file <- paste0(prediction_path, method, "-", gsub(" ", "", reg), "-", season, "-loso-predictions.rds")
        system(paste("cp", orig_file, new_file))
    }
    
}
