## KDE estimation
## 10/7/2016 - Nicholas Reich - file created
## 1/10/2017 - Evan Ray - moved to awes repository, minor modifications
## 1/11/2017 - Nicholas Reich - modified to create loso fits/predictions
## 1/16/2017 - Nicholas Reich - modified for standardized LOSO prediction function

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
foreach(region, analysis_time_season)
get_log_scores_via_direct_simulation(
    first_test_season = "2011/2012",
    all_seasons_left_out = paste0(1997:2010, "/", 1998:2011),
    analysis_time_season = analysis_time_season,
    first_analysis_time_season_week = 10, # == week 40 of year
    last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    region = region,
    prediction_target_var = "weighted_ili",
    incidence_bins = data.frame(
        lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
        upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1)),
    n_trajectory_sims = 100000,
    model_name = "kde",
    prediction_save_path = "inst/estimation/loso-predictions"
)

source("inst/estimation/kde/check-kde-predictions.R")
