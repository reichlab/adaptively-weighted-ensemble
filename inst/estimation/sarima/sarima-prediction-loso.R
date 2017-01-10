library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(MMWRweek)
library(awes)

### Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
region <- args[1]
analysis_time_season <- args[2]

## Parameters used in simulating trajectories via kcde
simulate_trajectories_sarima_params <- list(
  fits_filepath = "inst/estimation/sarima/fits",
  prediction_target_var = "weighted_ili",
  seasonal_difference = TRUE,
  transformation = "log"
)

get_log_scores_via_trajectory_simulation(
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
  simulate_trajectories_function = sample_predictive_trajectories_arima_wrapper,
  simulate_trajectories_params = simulate_trajectories_sarima_params,
  model_name = "sarima",
  prediction_save_path = "inst/estimation/loso-predictions"
)
