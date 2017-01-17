#' Simulate predictive distributions from KCDE/copulas model; suitable for use
#' as the \code{simulate_trajectories_function} argument to
#' \code{get_log_scores_via_trajectory_simulation}.
#' 
#' @param n_sims number of trajectories to simulate
#' @param max_prediction_horizon how many steps ahead to simulate
#' @param data data set
#' @param region region
#' @param analysis_time_season season in which we're predicting
#' @param analysis_time_season_week week of the season in which we're making our
#'   predictions, using all data up to the analysis time to make predictions for
#'   later time points
#' @param params other parameters.  A list with the following entries:
#'   * n_kcde_sims = number of simulations from kcde predictive distributions;
#'       this is different from number of predictive trajectories
#'   * copula_save_path = path to directory where copula fits are stored
#'   * estimation_results_path = path to directory where kcde fits are stored
#'   * max_lag = character giving maximum lag value used in kcde fits
#'   * seasonality = logical giving whether a seasonal kernel was used in kcde
#'   * bw_parameterization = string giving parameterization of bandwidth in kcde
#'       fit
#'   * last_analysis_time_season_week = 41
#' 
#' @return an n_sims by h matrix with simulated values
#' 
#' @export
simulate_trajectories_kcde <- function(
  n_sims,
  max_prediction_horizon,
  data,
  region,
  analysis_time_season,
  analysis_time_season_week,
  params
) {
  require(forecast) # for BoxCox transformations
  
  if(as.integer(substr(analysis_time_season, 1, 4)) >=
      as.integer(substr(params$first_test_season, 1, 4))) {
    ## if we're doing prediction during the test period,
    ## get fit obtained using all training data
    analysis_time_season <- "none"
  }
  
  ## load saved box-cox transformation parameters and do box-cox transformation
  lambda_case_descriptor <- case_descriptor <- paste0(
    region,
    "-prediction_horizon_", 1L,
    "-max_lag_", params$max_lag,
    "-seasonality_", params$seasonality,
    "-season_left_out_", gsub("/", "-", analysis_time_season)
  )
  
  lambda <- readRDS(
    file = file.path(params$estimation_results_path,
      paste0("box_cox_lambda-",
        lambda_case_descriptor,
        ".rds")
    )
  )
  
  data$box_cox_trans_weighted_ili <- BoxCox(data$weighted_ili, lambda)
  
  ## Load copula fits
  copula_fits_case_descriptor <- paste0(
    region,
    "-max_lag_", params$max_lag,
    "-seasonality_", params$seasonality,
    "-season_left_out_", gsub("/", "-", analysis_time_season)
  )
  
  copula_fits_file_name <- paste0(
    "kcde-copula-fits-",
    copula_fits_case_descriptor,
    ".rds"
  )
  
  copula_fits <- readRDS(file = file.path(params$copula_save_path, copula_fits_file_name))
  analysis_time_season_week_by_copula_fit <- unlist(lapply(copula_fits,
    function(copula_fit) {copula_fit$analysis_time_season_week}))
  
  kcde_fits_by_prediction_horizon <- lapply(seq_len(35),
    function(prediction_horizon) {
      case_descriptor <- paste0(
        region,
        "-prediction_horizon_", prediction_horizon,
        "-max_lag_", params$max_lag,
        "-seasonality_", params$seasonality,
        "-season_left_out_", gsub("/", "-", analysis_time_season)
      )
      
      kcde_fit_file_path <- file.path(params$estimation_results_path,
        paste0("kcde_fit-", case_descriptor, ".rds"))
      readRDS(kcde_fit_file_path)
    }
  )
  
  
  ## get the right copula for analysis_time_season_week
  predictive_copula_ind <- which(analysis_time_season_week_by_copula_fit == analysis_time_season_week)
  if(length(predictive_copula_ind > 0)) {
    copula_fit <- copula_fits[[predictive_copula_ind]]$copula_fit
    predictive_copula <- copula_fit@copula
  }
  
  sim_sequences <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
  na_rows <- seq_len(n_sims)
  for(sim_ind in na_rows) {
    ## simulate sequence from copula
    ## in the last few weeks of the year, we only needed a copula
    ## for trajectories that carry us through the end of the CDC
    ## season.  However, we may need predictions at individual weeks
    ## beyond that horizon.  In that case, pad with runif(), which
    ## is equivalent to assuming independence across horizons
    
    ## trajectories from copula
    if(length(predictive_copula_ind > 0)) {
      temp <- rCopula(1, predictive_copula)[1, ]
    } else {
      if(analysis_time_season_week != params$last_analysis_time_season_week) {
        stop(paste0(
          "Error: missing copula fit for analysis time season week ",
          analysis_time_season_week))
      }
      ## We don't have a copula fit because it's the last analysis time season week
      temp <- c()
    }
    
    ## if necessary, pad
    num_missing <- ncol(sim_sequences) - length(temp)
    if(num_missing > 0) {
      temp <- c(temp, runif(num_missing, min = 0, max = 1))
    }
    
    sim_sequences[sim_ind, ] <- temp
  }
  
  ### get quantiles from marginal predictive distributions corresponding to
  ### values simulated from copula
  trajectory_samples <- matrix(NA, nrow = n_sims, ncol = max_prediction_horizon)
  for(prediction_horizon in seq_len(max_prediction_horizon)) {
    trajectory_samples[, prediction_horizon] <-
      kcde_predict(
        p = sim_sequences[, prediction_horizon],
        n = params$n_kcde_sims,
        kcde_fit = kcde_fits_by_prediction_horizon[[prediction_horizon]],
        prediction_data = data,
        leading_rows_to_drop = 0L,
        trailing_rows_to_drop = 0L,
        additional_training_rows_to_drop = NULL,
        prediction_type = "quantile",
        log = TRUE
      )
  }
  
  return(trajectory_samples)
}
