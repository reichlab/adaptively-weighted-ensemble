## utility functions that might be useful for multiple prediction methods

#' Utility function to compute onset week based on a trajectory of incidence values
#' 
#' @param incidence_trajectory a numeric vector with incidence for each time
#'   point in a season
#' @param baseline the threshold that incidence must cross to count as an onset
#' @param onset_length number of consecutive time points that incidence must
#'   exceed the baseline threshold in order to count as the season onset
#' @param first_season_week number of weeks in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' @param weeks_in_first_season_year How many MMWR weeks are in the first year
#'   of the season?  For example, in the 2000/2001 season, the first year is
#'   2000.  There were 52 MMWR weeks in 2000.
#' 
#' @return the smallest index i such that every entry of
#'   incidence_trajectory[seq(from = i, length = onset_length)]
#'   is >= baseline, if such an index exists.
#'   Otherwise, the character vector "none"
#' 
#' @export
get_onset_week <- function(incidence_trajectory,
  baseline,
  onset_length,
  first_season_week = 31,
  weeks_in_first_season_year) {
  
  exceeded_threshold <- sapply(
    seq_len(length(incidence_trajectory) - onset_length),
    function(start_ind) {
      above_baseline <- incidence_trajectory[seq(from = start_ind, length = onset_length)] >= baseline
      length(above_baseline)>0 &&
        all(above_baseline) &&
        !all(is.na(incidence_trajectory))
    }
  )
  
  if(any(exceeded_threshold, na.rm = TRUE)) {
    season_week <- min(which(exceeded_threshold))
    
    return(season_week)
  } else {
    return("none")
  }
}

#' Convert season week to year week
#' 
#' @param season_week vector of indices of weeks in season
#' @param first_season_week number of week in year corresponding to the first
#'   week in the season.  For example, our code takes this value to be 31:
#'   a new influenza season starts on the 31st week of each year.
#' @param weeks_in_first_season_year How many MMWR weeks are in the first year
#'   of the season?  For example, in the 2000/2001 season, the first year is
#'   2000.  There were 52 MMWR weeks in 2000.
#' 
#' @return vector of the same length as season_week with the week of the year
#'   that each observation falls on
#' 
#' @export
season_week_to_year_week <- function(
  season_week,
  first_season_week = 31,
  weeks_in_first_season_year) {
  
  year_week <- season_week
  
  ## For competition first bin is week 40
  year_week[season_week < 10] <- 40
  year_week[season_week >= 10] <- season_week + first_season_week - 1
  year_week[year_week > weeks_in_first_season_year] <-
    year_week[year_week > weeks_in_first_season_year] -
    weeks_in_first_season_year
  
  return(year_week)
}

#' Convert year week to season week
#' 
#' @param year_week vector of indices of weeks in year
#' @param year either a single (four digit) year or a vector of years with the
#'   same length as year_week
#' 
#' @return vector of the same length as year_week with the week of the season
#'   that each observation falls on
#' 
#' @export
year_week_to_season_week <- function(
  year_week,
  year) {
  season_week <- ifelse(
    year_week <= 30,
    year_week + MMWRweek::MMWRweek(MMWRweek:::start_date(year) - 1)$MMWRweek - 30,
    year_week - 30
  )
  
  return(season_week)
}


#' Compute season onset, peak week, and peak incidence
#' 
#' @param data a data frame containing at minimum columns named season, 
#'   season_week and a column with some sort of incidence measure
#' @param season the season to look at
#' @param first_CDC_season_week the first week of the season to use for
#'   calculating onset and peak
#' @param last_CDC_season_week the last week of the season to use for
#'   calculating onset and peak
#' @param onset_baseline numeric baseline value for determining season onset
#' @param incidence_var a character string naming the variable in the data
#'   argument containing a measure of incidence, or an integer index
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#'   
#' @return a list with four entries:
#'   1) observed_onset_week, either an integer between first_CDC_season_week
#'     and last_CDC_season_week (inclusive), or "none"
#'   2) observed_peak_week, an integer between first_CDC_season_week and
#'     last_CDC_season_week (inclusive)
#'   3) observed_peak_inc, a numeric with the maximum value of the specified
#'     incidence measure between first_CDC_season_week and last_CDC_season_week
#'   4) observed_peak_inc_bin, character name of incidence bin for peak incidence
#' 
#' @export
get_observed_seasonal_quantities <- function(
  data,
  season,
  first_CDC_season_week = 10,
  last_CDC_season_week = 42,
  onset_baseline,
  incidence_var,
  incidence_bins,
  incidence_bin_names
) {
  first_season_ind <- min(which(data$season == season))
  last_season_ind <- max(which(data$season == season))
  
  obs_inc_in_season_leading_trailing_nas <-
    data[seq(from = first_season_ind, to = last_season_ind),
      incidence_var]
  
  ## pad so that we start at season week 1
  if(data$season_week[first_season_ind] != 1) {
    obs_inc_in_season_leading_trailing_nas <- c(
      rep(NA, data$season_week[first_season_ind] - 1),
      obs_inc_in_season_leading_trailing_nas)
  }
  
  ## set values before first analysis time season week or after last
  ## analysis time season week to NA
  ## these are outside of the bounds of the season the CDC wants to look at
  obs_inc_in_season_leading_trailing_nas[
    seq_len(first_CDC_season_week - 1)] <- NA
  if(length(obs_inc_in_season_leading_trailing_nas) >
      last_CDC_season_week) {
    obs_inc_in_season_leading_trailing_nas[
      seq(from = last_CDC_season_week + 1,
        to = length(obs_inc_in_season_leading_trailing_nas))] <- NA
  }
  
  observed_peak_inc <- max(
    obs_inc_in_season_leading_trailing_nas,
    na.rm = TRUE)
  observed_peak_inc_bin <- get_inc_bin(observed_peak_inc, return_character = TRUE)
  
  ## peak week timing is based on rounded values
  round_to_.1 <- function(inc_val) {
    if(is.na(inc_val)) {
      return(inc_val)
    } else {
      floor_val <- floor(inc_val * 10) / 10
      if(inc_val >= floor_val + 0.05) {
        return(floor_val + 0.1)
      } else {
        return(floor_val)
      }
    }
  }
  
  rounded_observed_peak_inc <- round_to_.1(observed_peak_inc)
  rounded_obs_inc_in_season <- sapply(obs_inc_in_season_leading_trailing_nas,
    round_to_.1
  )

  observed_peak_week <-
    which(rounded_obs_inc_in_season == as.numeric(rounded_observed_peak_inc))
  
  observed_onset_week <- get_onset_week(
    incidence_trajectory = rounded_obs_inc_in_season,
#    incidence_trajectory = obs_inc_in_season_leading_trailing_nas, # used in stable method
    baseline = onset_baseline,
    onset_length = 3L,
    first_season_week = 31,
    weeks_in_first_season_year = weeks_in_first_season_year
  )
  
  return(list(observed_onset_week = observed_onset_week,
    observed_peak_week = observed_peak_week,
    observed_peak_inc = observed_peak_inc,
    observed_peak_inc_bin = observed_peak_inc_bin
  ))
}


#' Calculate "log scores" for the purpose of the competition -- log[sum_i(p_i)] where p_i is the model's
#' probability of bin i and i runs over some bins adjacent to the bin where the observed quantity was.
#' 
#' @param bin_log_probs named numeric vector with log probability of each bin;
#'   names identify the bins
#' @param observed_bin character vector with name(s) of the observed bin(s)
#'   (Note that peak can occur in multiple bins)
#' @param prediction_target 
#' 
#' @return log score for the given observation
#' 
#' @export
compute_competition_log_score <- function(bin_log_probs,
  observed_bin,
  prediction_target) {
  ## validate probabilities sum to 1 and if not, force them to, with warning
  if(sum(exp(bin_log_probs)) != 1) {
    warning(paste(prediction_target, "probabilities do not sum to 1. automatically adjusting."))   
    bin_probs <- exp(bin_log_probs)
    bin_log_probs <- log(bin_probs/sum(bin_probs))
  }
    
  ## validate bin names match expected for prediction_target and
  ## observed_bin has appropriate length
  if(identical(prediction_target, "onset_week")) {
    expected_bin_names_52 <- c(as.character(10:42), "none")
    expected_bin_names_53 <- c(as.character(10:43), "none")
    
    if(!identical(length(observed_bin), 1L) || !identical(typeof(observed_bin), "character")) {
      stop("For prediction target onset_week, observed_bin must be a character vector of length 1")
    }
    
    if(identical(sort(expected_bin_names_52), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_52
    } else if(identical(sort(expected_bin_names_53), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_53
    } else {
      stop("invalid names for the vector bin_log_probs")
    }
  } else if(identical(prediction_target, "peak_week")) {
    expected_bin_names_52 <- as.character(10:42)
    expected_bin_names_53 <- as.character(10:43)
    
    if(length(observed_bin) == 0 || !identical(typeof(observed_bin), "character")) {
      stop("For prediction target onset_week, observed_bin must be a character vector of length > 0")
    }
    
    if(identical(sort(expected_bin_names_52), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_52
    } else if(identical(sort(expected_bin_names_53), sort(names(bin_log_probs)))) {
      expected_bin_names <- expected_bin_names_53
    } else {
      stop("invalid names for the vector bin_log_probs")
    }
  } else if(prediction_target %in% c("peak_inc", paste0("ph", 1:4, "_inc"))) {
    expected_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    if(!identical(length(observed_bin), 1L) || !identical(typeof(observed_bin), "character")) {
      stop("For given prediction target, observed_bin must be a character vector of length 1")
    }
    
    if(!identical(sort(expected_bin_names), sort(names(bin_log_probs)))) {
      stop("invalid names for the vector bin_log_probs")
    }
  } else {
    stop("Invalid value for prediction_target: must be 'onset_week', ")
  }
  
  ## validate observed_bin is one of the expected_bin_names
  if(!(all(observed_bin %in% expected_bin_names))) {
    stop(paste0(
      "observed_bin must be one of (",
      paste(expected_bin_names, collapes = ", "),
      ")"
    ))
  }
  
  ## get bins to sum over
  obs_bin_inds <- sapply(observed_bin, function(bin_name) {
    which(expected_bin_names == bin_name)
  })
  if(identical(prediction_target, "onset_week")) {
    if(identical(observed_bin, "none")) {
      bins_to_sum <- obs_bin_inds
    } else {
      bins_to_sum <- obs_bin_inds + seq(from = -1, to = 1, by = 1)
      bins_to_sum <- bins_to_sum[
        bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names) - 1]
    }
  } else if(identical(prediction_target, "peak_week")) {
    bins_to_sum <- unique(as.vector(sapply(obs_bin_inds, function(bin_ind) {
      bin_ind + seq(from = -1, to = 1, by = 1)
    })))
    bins_to_sum <- bins_to_sum[
      bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names)]
  } else if(prediction_target %in% c("peak_inc", paste0("ph", 1:4, "_inc"))) {
    bins_to_sum <- obs_bin_inds + seq(from = -5, to = 5)
    bins_to_sum <- bins_to_sum[
      bins_to_sum >= 1 & bins_to_sum <= length(expected_bin_names)]
  }
  
  ## Do summation
  ## Futz around with bin names because order of expected bin names may not
  ## match order of bin_log_probs
  bin_names_to_sum <- expected_bin_names[bins_to_sum]
  log_prob <- logspace_sum(bin_log_probs[bin_names_to_sum])
  
  ## They truncate at -10
  log_prob <- max(-10, log_prob)
  
  ## return
  return(log_prob)
}

#' Get the onset baseline for a combination of region and season
#' 
#' @param region a string, either "X", "Region k", or "Regionk" where
#'   k \in {1, ..., 10}
#' @param season a string, in the format "2015/2016"
#' 
#' @return baseline value for determining season onset
#' 
#' @export
get_onset_baseline <- function(region, season = "2015/2016") {
  ## pick baseline
  baselines <- read.csv(file="data-raw/cdc-baselines.csv")
  
  ## assumes region is either "X" or "Region k" format
  reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
  idx <- which(baselines$region==reg_string & baselines$season==season)
  reg_baseline <- baselines[idx, "baseline"]
  
  return(reg_baseline)
}

#' Helper function to standardize creation of empty prediction dataset
#'
#' @param data the data to be used
#' @param model_name the name of the model 
#' @param incidence_bin_names character vector of names for incidence bins
#' @param first_analysis_time_season_week the first week during each season in
#'   which we will form a prediction
#' @param last_analysis_time_season_week the last week during each season in
#'   which we will form a prediction
#'
#' @return the prediction_df filled with NAs
#' 
#' @export
make_predictions_dataframe <- function(data,
  model_name,
  incidence_bin_names,
  first_analysis_time_season_week = 10,
  last_analysis_time_season_week = 41) {
  
  ## allocate more than enough space up front,
  ## delete extra later
  na_vec <- rep(NA,
    (last_analysis_time_season_week - first_analysis_time_season_week + 1)
  )
  
  onset_week_bins <- c(as.character(10:42), "none")
  peak_week_bins <- as.character(10:42)
  
  predictions_df <- cbind(
    data.frame(
      model = model_name,
      analysis_time_season = na_vec,
      analysis_time_season_week = na_vec,
      prediction_week_ph_1 = na_vec,
      prediction_week_ph_2 = na_vec,
      prediction_week_ph_3 = na_vec,
      prediction_week_ph_4 = na_vec,
      onset_log_score = na_vec,
      peak_week_log_score = na_vec,
      peak_inc_log_score = na_vec,
      ph_1_inc_log_score = na_vec,
      ph_2_inc_log_score = na_vec,
      ph_3_inc_log_score = na_vec,
      ph_4_inc_log_score = na_vec,
      onset_competition_log_score = na_vec,
      peak_week_competition_log_score = na_vec,
      peak_inc_competition_log_score = na_vec,
      ph_1_inc_competition_log_score = na_vec,
      ph_2_inc_competition_log_score = na_vec,
      ph_3_inc_competition_log_score = na_vec,
      ph_4_inc_competition_log_score = na_vec,
      stringsAsFactors = FALSE), # don't want model to be a factor with 1 level
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(onset_week_bins))) %>%
      `colnames<-`(paste0("onset_bin_", onset_week_bins, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(peak_week_bins))) %>%
      `colnames<-`(paste0("peak_week_bin_", peak_week_bins, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("peak_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_1_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_2_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_3_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_4_inc_bin_", incidence_bin_names, "_log_prob"))
  )
  
  return(predictions_df)
}

#' Get log scores and full predictive distributions for each prediction target
#' 
#' for a given season using a predictive method that works by simulating
#' trajectories of incidence in each remaining week of the season.  Results are
#' stored in a data frame, saved in a .rds file with a name like
#' "model_name-region-season-loso-predictions.rds"
#' Results have columns indicating the analysis time season and season week,
#' model name, log scores for each prediction target, the "log score" used
#' in the competition (adding probabilities from adjacent bins) for each
#' prediction target, as well as the log of the probability assigned to each
#' bin.
#' 
#' @param all_seasons_left_out character vector of seasons that were in the
#'   training period, in the format "2000/2001"
#' @param analysis_time_season character vector of length 1 specifying the
#'   season to obtain predictions for, in the format "2000/2001"
#' @param first_analysis_time_season_week integer specifying the first week of
#'   the season in which to make predictions, using all data up to and
#'   including that week to make predictions for each following week in the
#'   season
#' @param last_analysis_time_season_week integer specifying the last week of
#'   the season in which to make predictions, using all data up to and including
#'   that week to make predictions for each following week in the season
#' @param region string specifying the region to use, in the format "X" or
#'   "Region k" where k \in {1, ..., 10}
#' @param prediction_target_var string specifying the name of the variable in
#'   data for which we want to make predictions
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#' @param n_trajectory_sims integer number of trajectories to simulate
#' @param simulate_trajectories_function a function to call to simulate
#'   incidence trajectories.  It will be called with the following arguments:
#'     * n_sims = number of trajectories to simulate
#'     * max_prediction_horizon = number of following weeks to simulate
#'     * data = all available data to use in doing simulation, up to and
#'         including the analysis time
#'     * region = region
#'     * analysis_time_season = analysis_time_season
#'     * analysis_time_season_week = week of the season at which we are making
#'         the predictions
#'     * params = simulate_trajectories_params; additional user-provided
#'         parameters
#' @param simulate_trajectories_params optional additional parameters to pass
#'   to simulate_trajectories_function
#' @param model_name name of model, stored in the results data frame
#' @param prediction_save_path path to directory where results will be saved
#' 
#' @return none
#' 
#' @export
get_log_scores_via_trajectory_simulation <- function(
  all_seasons_left_out,
  analysis_time_season,
  first_analysis_time_season_week = 10, # == week 40 of year
  last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
  region,
  prediction_target_var,
  incidence_bins,
  incidence_bin_names,
  n_trajectory_sims,
  simulate_trajectories_function,
  simulate_trajectories_params,
  model_name,
  prediction_save_path
) {
  ## Load data.  The only reason to do this here is to know what the dimensions
  ## of the results data frame should be.
  data <- read.csv("data-raw/allflu-cleaned.csv", stringsAsFactors = FALSE)
  
  data$time <- as.POSIXct(data$time)
  
  ## subset data to be only the region of interest
  data <- data[data$region == region,]
  
  ## data frame to describe predictions
  ## allocate more than enough space up front,
  ## delete extra later
  predictions_df <- make_predictions_dataframe(
    data = data,
    model_name = model_name,
    incidence_bin_names = incidence_bin_names,
    first_analysis_time_season_week = 10,
    last_analysis_time_season_week = 41)
  
  results_save_row <- 1L
  
  ## make predictions for each prediction target in the left-out season
  ## for each possible "last observed" week, starting with the last week of the previous season
  
  ## get observed quantities related to overall season, for computing log scores
  observed_seasonal_quantities <- get_observed_seasonal_quantities(
    data = data,
    season = analysis_time_season,
    first_CDC_season_week = first_analysis_time_season_week,
    last_CDC_season_week = last_analysis_time_season_week + 1,
    onset_baseline = 
      get_onset_baseline(region = region, season = analysis_time_season),
    incidence_var = prediction_target_var,
    incidence_bins = incidence_bins,
    incidence_bin_names = incidence_bin_names
  )
  
  last_analysis_time_season_week_in_data <- max(data$season_week[data$season == analysis_time_season])
  for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = min(last_analysis_time_season_week, last_analysis_time_season_week_in_data - 1))) {
    cat(paste0("analysis season week = ", analysis_time_season_week, "\n"))
    
    analysis_time_ind <- which(data$season == analysis_time_season &
        data$season_week == analysis_time_season_week)
    
    ## keep track of if we made any predictions with this combination of season and week
    made_predictions <- FALSE
    
    ## simulate n_trajectory_sims trajectories
    max_prediction_horizon <- max(4L,
      last_analysis_time_season_week + 1 - analysis_time_season_week)
    
    trajectory_samples <- simulate_trajectories_function(
      n_sims = n_trajectory_sims,
      max_prediction_horizon = max_prediction_horizon,
      data = data[seq_len(analysis_time_ind), , drop = FALSE],
      region = region,
      analysis_time_season = analysis_time_season,
      analysis_time_season_week = analysis_time_season_week,
      params = simulate_trajectories_params
    )
    
    ## Round to nearest 0.1 -- they do this in competition
    ## I'm not using R's round() function because I want
    ## 0.05 -> 0.1
    ## (...?  We should see how they round?  Maybe this never comes up?
    ## I have a feeling I'm fretting over a probability 0 event.)
    trajectory_samples <- apply(
      trajectory_samples,
      c(1, 2),
      function(inc_val) {
        if(is.na(inc_val)) {
          return(inc_val)
        } else {
          floor_val <- floor(inc_val * 10) / 10
          if(inc_val >= floor_val + 0.05) {
            return(floor_val + 0.1)
          } else {
            return(floor_val)
          }
        }
      }
    )
    
    ## get indices in trajectory samples with NA values that affect
    ## estimation of seasonal quantities
    sample_inds_with_na <- apply(
      trajectory_samples[,
        seq_len(last_analysis_time_season_week + 1 - analysis_time_season_week),
        drop = FALSE],
      1,
      function(x) any(is.na(x)))
    
    ## Predictions for things about the whole season
    if(!all(sample_inds_with_na)) {
      made_predictions <- TRUE
      
      ## subset to sampled trajectories that are usable/do not have NAs
      subset_trajectory_samples <- trajectory_samples[!sample_inds_with_na, ]
      
      ## Augment trajectory samples with previous observed incidence values
      ## This is where we should be adding in something to account for backfill.
      first_season_obs_ind <- min(which(data$season == analysis_time_season))
      subset_trajectory_samples <- cbind(
        matrix(
          rep(data[seq(from = first_season_obs_ind, to = analysis_time_ind), prediction_target_var], each = n_trajectory_sims),
          nrow = nrow(subset_trajectory_samples)
        ),
        subset_trajectory_samples
      )
      
      ## If first observation for the season was not at season week 1,
      ## augment with leading NAs
      first_season_obs_week <- data$season_week[first_season_obs_ind]
      if(first_season_obs_week != 1) {
        subset_trajectory_samples <- cbind(
          matrix(NA, nrow = nrow(subset_trajectory_samples), ncol = first_season_obs_week - 1),
          subset_trajectory_samples
        )
      }
      
      ## values before the first analysis time week are NA so that
      ## onset and peak calculations only look at data within the CDC's definition
      ## of the flu season for purposes of the competition
      subset_trajectory_samples[,
        seq(from = 1, to = first_analysis_time_season_week - 1)] <- NA
      
      ## Convert to binned values
      binned_subset_trajectory_samples <- 
        get_inc_bin(subset_trajectory_samples,
          return_character = FALSE)
      
      ## Get onset week for each simulated trajectory
      onset_week_by_sim_ind <-
        apply(binned_subset_trajectory_samples, 1, function(trajectory) {
          get_onset_week(
            incidence_trajectory = trajectory,
            baseline = 
              get_onset_baseline(region = region, season = analysis_time_season),
            onset_length = 3L,
            first_season_week = 31,
            weeks_in_first_season_year = weeks_in_first_season_year
          )
        })
      
      ## Get peak incidence for each simulated trajectory
      peak_inc_bin_by_sim_ind <-
        apply(binned_subset_trajectory_samples, 1, function(trajectory) {
          max(trajectory, na.rm = TRUE)
        })
      
      ## get peak week by sim ind
      ## note that some sim inds may have more than 1 peak week...
      peak_weeks_by_sim_ind <- unlist(lapply(
        seq_len(nrow(binned_subset_trajectory_samples)),
        function(sim_ind) {
          bin_val <- peak_inc_bin_by_sim_ind[sim_ind]
          peak_season_weeks <- which(
            binned_subset_trajectory_samples[sim_ind, ] == bin_val)
          return(peak_season_weeks)
        }
      ))
      
      ## Get log scores
      onset_week_bins <- c(as.character(10:42), "none")
      onset_bin_log_probs <- log(sapply(
        onset_week_bins,
        function(bin_name) {
          sum(onset_week_by_sim_ind == bin_name)
        })) -
        log(length(onset_week_by_sim_ind))
      predictions_df[results_save_row, paste0("onset_bin_", onset_week_bins, "_log_prob")] <-
        onset_bin_log_probs
      predictions_df[results_save_row, "onset_log_score"] <-
        onset_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_onset_week)]
      predictions_df[results_save_row, "onset_competition_log_score"] <-
        compute_competition_log_score(onset_bin_log_probs,
          as.character(observed_seasonal_quantities$observed_onset_week),
          "onset_week")
      
      peak_week_bins <- as.character(10:42)
      peak_week_bin_log_probs <- log(sapply(
        peak_week_bins,
        function(bin_name) {
          sum(peak_weeks_by_sim_ind == bin_name)
        })) -
        log(length(peak_weeks_by_sim_ind))
      names(peak_week_bin_log_probs) <- as.character(peak_week_bins)
      predictions_df[results_save_row, paste0("peak_week_bin_", peak_week_bins, "_log_prob")] <-
        peak_week_bin_log_probs
      predictions_df[results_save_row, "peak_week_log_score"] <-
        logspace_sum(peak_week_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_peak_week)])
      predictions_df[results_save_row, "peak_week_competition_log_score"] <-
        compute_competition_log_score(peak_week_bin_log_probs,
          as.character(observed_seasonal_quantities$observed_peak_week),
          "peak_week")
      
      peak_inc_bin_log_probs <- log(sapply(
        incidence_bin_names,
        function(bin_name) {
          sum(peak_inc_bin_by_sim_ind == as.numeric(bin_name))
        })) -
        log(length(peak_inc_bin_by_sim_ind))
      predictions_df[results_save_row, paste0("peak_inc_bin_", incidence_bin_names, "_log_prob")] <-
        peak_inc_bin_log_probs
      predictions_df[results_save_row, "peak_inc_log_score"] <-
        peak_inc_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_peak_inc_bin)]
      predictions_df[results_save_row, "peak_inc_competition_log_score"] <-
        compute_competition_log_score(peak_inc_bin_log_probs,
          observed_seasonal_quantities$observed_peak_inc_bin,
          "peak_inc")
    }
    
    ## Predictions for incidence in an individual week at prediction horizon ph = 1, ..., 4
    for(ph in 1:4) {
      sample_inds_with_na <- is.na(trajectory_samples[, ph])
      
      ## get observed value/bin
      observed_ph_inc <-
        data[analysis_time_ind + ph, prediction_target_var]
      observed_ph_inc_bin <- get_inc_bin(observed_ph_inc, return_character = TRUE)
      
      if(!all(sample_inds_with_na) && !is.na(observed_ph_inc)) {
        made_predictions <- TRUE
        
        ## get sampled incidence values at prediction horizon that are usable/not NAs
        ph_inc_by_sim_ind <- trajectory_samples[!sample_inds_with_na, ph]
        ph_inc_bin_by_sim_ind <- get_inc_bin(ph_inc_by_sim_ind, return_character = TRUE)
        
        ## get log score
        ph_inc_bin_log_probs <- log(sapply(
          incidence_bin_names,
          function(bin_name) {
            sum(ph_inc_bin_by_sim_ind == bin_name)
          })) -
          log(length(ph_inc_bin_by_sim_ind))
        predictions_df[results_save_row, paste0("ph_", ph, "_inc_bin_", incidence_bin_names, "_log_prob")] <-
          ph_inc_bin_log_probs
        predictions_df[results_save_row, paste0("ph_", ph, "_inc_log_score")] <-
          ph_inc_bin_log_probs[observed_ph_inc_bin]
        predictions_df[results_save_row,
          paste0("ph_", ph, "_inc_competition_log_score")] <-
          compute_competition_log_score(ph_inc_bin_log_probs,
            observed_ph_inc_bin,
            paste0("ph", ph, "_inc"))
      }
    } # ph loop
    
    if(made_predictions) {
      predictions_df[results_save_row, "analysis_time_season"] <- analysis_time_season
      predictions_df[results_save_row, "analysis_time_season_week"] <- analysis_time_season_week
      predictions_df[results_save_row, "prediction_week_ph_1"] <- analysis_time_season_week + 1
      predictions_df[results_save_row, "prediction_week_ph_2"] <- analysis_time_season_week + 2
      predictions_df[results_save_row, "prediction_week_ph_3"] <- analysis_time_season_week + 3
      predictions_df[results_save_row, "prediction_week_ph_4"] <- analysis_time_season_week + 4
      
      results_save_row <- results_save_row + 1
    }
  } # analysis_time_season_week
  
  ## if there are extra rows in the predictions_df, delete them
  if(results_save_row <= nrow(predictions_df)) {
    predictions_df <- predictions_df[
      -seq(from = results_save_row, to = nrow(predictions_df)),
      ,
      drop = FALSE
      ]
  }
  
  region_str <- ifelse(identical(region, "X"), "National", gsub(" ", "", region))
  season_str <- gsub("/", "-", analysis_time_season)
  saveRDS(predictions_df,
    file = paste0(prediction_save_path,
      "/",
      model_name, "-", region_str, "-", season_str, "-loso-predictions.rds"))
}


#' return the bin name for a given incidence
#' 
#' @param inc numeric incidence level
#' @param return_character logical: if true, return type is character (bin name)
#'   if false, return type is numeric representation of bin
#' 
#' @return vector giving the bin name of the input incidence.
#' 
#' @details assumes max inc bin is 13 and bins are 0.1 in size. 
#' 
#' @export
get_inc_bin <- function(inc,
    return_character = TRUE) {
  inc <- round(inc, 1)
  bin_numeric <- ifelse(inc < 13,
    floor(inc*10)/10, ## floors to 1st decimal place
    13)
  if(return_character) {
    return(as.character(bin_numeric))
  } else {
    return(bin_numeric)
  }
}


#' Calcluation of median value from binned probability distribution
#' 
#' @param probs vector of named probabilities
#' 
#' @return a numeric value
#' 
#' @export
calc_median_from_binned_probs <- function(probs) {
    ## could do something more intelligent for "none" bin in onset - currently assuming it is all ordered
    cumprob <- cumsum(probs)
    median_idx <- min(which(cumprob>=0.5))
    as.numeric(names(probs)[median_idx])
}


#' Get the initial rng substream for an rstream object.  Should be called by all
#' prediction methods before doing any random number generation.
#' 
#' This function DOES NOT set RNG to use rstream with the returned object.
#' Because of strange behavior in the rstream package, this function can be
#' called only once in a given R session.
#' 
#' @param seed integer seed for rng; the default was randomly generated
#' 
#' @return object of class "rstream.mrg32k3a".  The object has been packed via
#'   rstream.packed.
#' 
#' @export
get_initial_rng_substream <- function(
  seed = 4114043) {
  require("rstream")
  
  set.seed(seed)
  rngstream <- new("rstream.mrg32k3a", seed = sample(1:100000, 6, rep = FALSE))
  
  ## pack rngstream object and return (invisibly) in case methods want to use
  rstream.packed(rngstream) <- TRUE
  return(rngstream)
}


#' Get the rng substream for an rstream object corresponding to the combination
#' of prediction method, region, and season left out.  Should be called by all
#' prediction methods before doing any random number generation.
#' 
#' Importantly, by default this function has the side effect of setting RNG to
#' use rstream with the returned object.  This behavior is determined by the
#' set_rng argument.  This means the caller doesn't have to worry about doing
#' anything unless (a) it wants to use more than 1 substream or (b) it is going
#' to parallelize or do any RNG in a different R session.  Because of strange
#' behavior in the rstream package, this function can be called at most once
#' without the rngstream argument in a given R session.
#' 
#' @param rngstream (optional) object of class "rstream.mrg32k3a" which will be
#'   advanced from its current state.
#' @param seed integer seed for rng; the default was randomly generated
#' @param method character string specifying prediction method
#'   currently one of "sarima", "kcde", or "kde"
#' @param region character string specifying region; format "National" or
#'   "Region k" for k \in {1, ..., 10}
#' @param season character string specifying season predicted, in format
#'   "1998/1999"
#' @param set_rng boolean should rng be set to use rstream with the returned
#'   rngstream object?
#' 
#' @return invisible object of class "rstream.mrg32k3a", but advanced to the
#' first substream reserved for the given combination of prediction method,
#' region, and season.  The object has been packed via rstream.packed.
#' 
#' @export
get_rng_substream <- function(
  rngstream,
  seed = 4114043,
  method,
  region,
  season,
  set_rng = TRUE) {
  require("rstream")
  
  ## Create a data frame with combinations of method, region, and season,
  ## number of substreams used for each such combination.
  ## We can add more methods later without causing any problems by appending
  ## new method names to the END of the "method" vector below.
  ## Adding new seasons or regions must be done by adding a new set of rows
  ## to the bottom of the substreams_used data frame (e.g. via bind_rows).
  substreams_used <- expand.grid(
    season = paste0(1997:2015, "/", 1998:2016),
    region = c("National", paste0("Region ", 1:10)),
    method = c("sarima", "kcde", "kde"),
    stringsAsFactors = FALSE
  )
  substreams_used$num_substreams <- 1
  ## if any future method uses more than 1 substream, set that here
  
  ## substream index for the specified method, region, and season
  if(identical(region, "X")) {
    region <- "National"
  }
  ind <- which(
    substreams_used$season == season &
    substreams_used$region == region &
    substreams_used$method == method)
  
  if(!identical(length(ind), 1L)) {
    stop("Invalid season, region, and/or method.")
  }
  
  ## Create Rstream object and advance past all substreams used by previous
  ## methods/regions/seasons
  if(missing(rngstream)) {
    set.seed(seed)
    rngstream <- new("rstream.mrg32k3a", seed = sample(1:100000, 6, rep = FALSE))
  } else {
    rstream.packed(rngstream) <- FALSE
  }
  
  advance_count <- sum(substreams_used$num_substreams[seq_len(ind - 1)])
  for(i in seq_len(advance_count)) {
    rstream.nextsubstream(rngstream)
  }
  
  ## set to use rstream package for RNG with rngstream object
  if(set_rng) {
    rstream.RNG(rngstream)
  }
  
  ## pack rngstream object and return (invisibly) in case methods want to use
  rstream.packed(rngstream) <- TRUE
  invisible(rngstream)
}


## interface to R's C API for logspace arithmetic

#' Calculate log(exp(logx) - exp(logy)) in a somewhat numerically stable way.
#' 
#' @param logx, logy log-scale numeric values to subtract
#' 
#' @return log(exp(logx) - exp(logy)), but more numerically stable
#' 
#' @export
logspace_sub <- function(logx, logy) {
  return(.Call("logspace_sub_C",
    as.numeric(logx),
    as.numeric(logy),
    PACKAGE = "awes"))
}

#' Calculate log(exp(logx) + exp(logy)) in a somewhat numerically stable way.
#' 
#' @param logx, logy log-scale numeric values to add
#' 
#' @return log(exp(logx) + exp(logy)), but more numerically stable
#' 
#' @export
logspace_add <- function(logx, logy) {
  return(.Call("logspace_add_C",
    as.numeric(logx),
    as.numeric(logy),
    PACKAGE = "awes"))
}

#' Calculate log(sum(exp(logx))) in a somewhat numerically stable way.
#' 
#' @param logx log-scale numeric vector of values to sum
#' 
#' @return log(sum(exp(logx))), but more numerically stable
#' 
#' @export
logspace_sum <- function(logx) {
  dim(logx) <- c(1, length(logx))
  return(logspace_sum_matrix_rows(logx))
}

#' Calculate logspace summation of matrix rows in a somewhat numerically stable
#' way.
#' 
#' @param logX log-scale numeric matrix of values to sum.
#' 
#' @return log(apply(exp(logX), 1, sum)), but more numerically stable
#' 
#' @export
logspace_sum_matrix_rows <- function(logX) {
  return(.Call("logspace_sum_matrix_rows_C",
    as.numeric(logX),
    as.integer(nrow(logX)),
    as.integer(ncol(logX)),
    PACKAGE = "awes"))
}

#' Calculate logspace difference of matrix rows in a somewhat numerically stable
#' way.
#' 
#' @param logX log-scale numeric matrix of values to subtract.  logX must have
#'   exactly 2 columns.
#' 
#' @return log(exp(logX)[, 1] - exp(logX)[, 2]), but more numerically stable
#' 
#' @export
logspace_sub_matrix_rows <- function(logX) {
  if(!is.matrix(logX) || !identical(ncol(logX), 2L))
    stop("logX must be a matrix with 2 columns")
  
  return(.Call("logspace_sub_matrix_rows_C",
    as.numeric(logX),
    as.integer(nrow(logX)),
    PACKAGE = "awes"))
}
