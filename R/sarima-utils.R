## utility functions for SARIMA fits

## This is directly taken from the forecast.Arima function from the forecast package,
## but I've forced bootstrap = TRUE and return the simulated trajectories which were
## not returned in the original function definition.
sample_predictive_trajectories_arima <- function (object,
    h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
    level = c(80, 95),
    fan = FALSE,
    xreg = NULL,
    lambda = object$lambda,
    npaths = 5000,
    ...) {
    sim <- matrix(NA, nrow = npaths, ncol = h)
    
    for (i in 1:npaths) {
        try(sim[i, ] <- simulate.Arima(object, nsim = h), silent = TRUE)
    }
    
    return(sim)
}

sample_predictive_trajectories_arima_wrapper <- function(
    n_sims,
    max_prediction_horizon,
    data,
    region,
    analysis_time_season,
    analysis_time_season_week,
    params
) {
    ## load SARIMA fit
    if(identical(region, "National") || identical(region, "X")) {
        reg_str <- "national"
    } else {
        space_ind <- grep(" ", strsplit(region, "")[[1]])
        reg_num <- as.integer(substr(region, space_ind + 1, nchar(region)))
        reg_str <- paste0("region", sprintf("%02d", reg_num))
    }
    
    if(identical(analysis_time_season, "2016/2017")) {
        ## if we're doing prediction for the actual competition,
        ## get fit obtained by leaving out an early season
        if(identical(reg_str, "region01") || identical(reg_str, "region08")) {
            fit_season_name <- "1999/2000" # for region 1, fit from 1998/1999 is nonstationary...
        } else {
            fit_season_name <- "1998/1999"
        }
    } else {
        fit_season_name <- analysis_time_season
    }
    fit_filepath <- file.path(
        params$fits_filepath,
        paste0(
            "sarima-",
            reg_str,
            "-fit-leave-out-",
            gsub("/", "-", fit_season_name),
            ".rds"))
    
    ## If no SARIMA fit, exit early by returning a matrix of NAs
    if(!file.exists(fit_filepath)) {
        return(matrix(NA, nrow = n_sims, ncol = max_prediction_horizon))
    }
    
    sarima_fit <- readRDS(file = fit_filepath)
    
    ## Update SARIMA fit object with seasonally differenced data up
    ## through analysis_time_ind
    if(identical(params$transformation, "log")) {
      undifferenced_new_target_data <- log(data[, params$prediction_target_var])
    } else {
      undifferenced_new_target_data <- data[, params$prediction_target_var]
    }
    if(params$seasonal_difference) {
      new_target_data <- ts(c(rep(NA, 52),
        undifferenced_new_target_data[seq(from = 53, to = nrow(data))] -
          undifferenced_new_target_data[seq(from = 1, to = nrow(data) - 52)]),
        frequency = 52)
    } else {
      new_target_data <- ts(undifferenced_new_target_data, frequency = 52)
    }

    ## deal with internal missing values in new_target_data that can result in 
    ## simulated trajectories of all NAs if the model has a moving average component
    ## 
    ## here we do this by linear interpolation
    ## 
    ## another solution would be to write a version of stats::filter that
    ## does the "right thing" with NAs if filter coefficients are 0
    ## and then use that function in forecast:::myarima.sim
    if(any(is.na(new_target_data)) && !is.na(tail(new_target_data, 1))) {
        ## drop leading NAs
        if(is.na(new_target_data[1])) {
            num_leading_nas <- rle(is.na(new_target_data))$lengths[1]
            new_target_data <- new_target_data[- seq_len(num_leading_nas)]
        }
        
        ## interpolate internal NAs
        while(any(is.na(new_target_data))) {
            na_rle <- rle(is.na(new_target_data))
            na_run_start_ind <- na_rle$lengths[1] + 1
            na_run_end_ind <- na_run_start_ind + na_rle$lengths[2] - 1
            new_target_data[na_run_start_ind:na_run_end_ind] <-
                approx(
                    x = c(na_run_start_ind - 1, na_run_end_ind + 1),
                    y = new_target_data[c(na_run_start_ind - 1, na_run_end_ind + 1)],
                    xout = na_run_start_ind:na_run_end_ind,
                    method = "linear"
                    )$y
        }
    }
    
    updated_sarima_fit <- Arima(
        new_target_data,
        model = sarima_fit)
    
    ## The update process via call to Arima can somehow result in 0/negative
    ## variance estimates.  Reset to original values
    updated_sarima_fit$sigma2 <- sarima_fit$sigma2
    updated_sarima_fit$var.coef <- sarima_fit$var.coef
    
    raw_trajectory_samples <- sample_predictive_trajectories_arima(
        updated_sarima_fit,
        h = max_prediction_horizon,
        npaths = n_sims)
    
    ## Sampled trajectories are of seasonally differenced log incidence
    ## Get to trajectories for originally observed incidence ("inc") by
    ## adding seasonal lag of incidence and exponentiating
    inc_trajectory_samples <- raw_trajectory_samples
    if(params$seasonal_difference) {
      for(prediction_horizon in seq_len(max_prediction_horizon)) {
        inc_trajectory_samples[, prediction_horizon] <-
          raw_trajectory_samples[, prediction_horizon] +
          undifferenced_new_target_data[nrow(data) + prediction_horizon - 52]
      }
    }
    if(identical(params$transformation, "log")) {
      inc_trajectory_samples <- exp(inc_trajectory_samples)
    }
    
    return(inc_trajectory_samples)
}


#' Estimate series of leave-one-season-out SARIMA models
#'
#' @param data regional dataset with structure like regionflu-cleaned
#' @param reg_num region number for estimation
#' @param first_test_season string indicating first test season
#' @param d order of first differencing
#' @param D order of seasonal differencing
#' @param seasonal_difference boolean; take a seasonal difference before passing
#'   to auto.arima?
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", or "none"
#' @param path path in which to save files
#'
#' @return NULL just saves files
#'
fit_region_sarima <- function(data, reg_num, first_test_season, d = NA, D = NA, seasonal_difference = TRUE, transformation = "none", path) {
    require(forecast)
    ## subset data to be only the region of interest
    if(reg_num == "X") {
        region_string <- "X"
    } else {
        region_string <- paste("Region", reg_num)
    }
    data <- data[data$region == region_string,]
    
    ## Subset data to do estimation using only data up through first_test_season
    ## remainder are held out for evaluating performance.
    first_ind_test_season <- min(which(data$season == first_test_season))
    data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]
    
    prediction_target_var <- "weighted_ili"
    
    if(identical(transformation, "log")) {
      undifferenced_target_data <- log(data[, prediction_target_var])
    } else {
      undifferenced_target_data <- data[, prediction_target_var]
    }
    
    if(seasonal_difference) {
      pred_target <- ts(c(rep(NA, 52),
        undifferenced_target_data[seq(from = 53, to = nrow(data))] -
          undifferenced_target_data[seq(from = 1, to = nrow(data) - 52)]),
        frequency = 52)
    } else {
      pred_target <- ts(undifferenced_target_data, frequency = 52)
    }
    
    for(season_left_out in c("none", unique(data$season))) {
        if(reg_num == "X") {
            filename <- paste0(
                path,
                "sarima-national",
                "-fit-leave-out-",
                gsub("/", "-", season_left_out),
                ".rds")
        } else {
            reg_num <- as.numeric(reg_num)
            filename <- paste0(
                path,
                "sarima-region",
                sprintf("%02d", reg_num),
                "-fit-leave-out-",
                gsub("/", "-", season_left_out),
                ".rds")
        }
        if(file.exists(filename)){
            message(paste("Region", reg_num, "leaving out", season_left_out, "already fit ::", Sys.time()))
            next
        }
        ## Estimate SARIMA model wihout data from given season, only if there would be
        ## anything to predict in that season if it were held out.
        if(identical(season_left_out, "none") || !all(is.na(pred_target[data$season == season_left_out]))) {
            message(paste("fitting Region", reg_num, "leaving out", season_left_out, "::", Sys.time()))
            pred_target_season_NAs <- pred_target
            pred_target_season_NAs[data$season == season_left_out] <- NA
            
            if(identical(transformation, "box-cox")) {
              lambda <- BoxCox.lambda(pred_target_season_NAs)
              
              ## get SARIMA fit
              sarima_fit_season_left_out <- auto.arima(pred_target_season_NAs,
                  d = d,
                  D = D,
                  stationary = TRUE,
                  lambda = lambda)
            } else {
              ## log transformation (already done manually) or no transformation
              sarima_fit_season_left_out <- auto.arima(pred_target_season_NAs,
                d = d,
                D = D,
                stationary = TRUE)
            }
            
            saveRDS(sarima_fit_season_left_out, file = filename)
        }
    }
}

#' Create regional predictions from a seasonal model
#'
#' @param data regionalflu-cleaned dataset
#' @param reg_num region number
#' @param first_test_season the first season reserved for testing data. seasons before this are predicted.
#' @param n_sim number of draws to simulate from each predictive distribution.
#' @param filepath path to folder with estimated fits
#'
#' @return NULL just saves a file
predict_region_sarima <- function(data, reg_num, first_test_season, n_sim, filepath) {
    require(forecast)
    
    ## load baselines for onset calculations
    baselines <- read.csv("data-raw/cdc-baselines.csv")
    baselines <- baselines[baselines$region == paste0("Region", reg_num),]
    
    ## Subset data to do estimation using only data up through first_test season
    ## subsequent years are held out for evaluating performance.
    first_ind_test_season <- min(which(data$season == first_test_season))
    data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]
    
    ## subset data to be only the region of interest
    region_string <- paste("Region", reg_num)
    data <- data[data$region == region_string,]
    
    prediction_target_var <- "weighted_ili"
    data$log_prediction_target <- log(data[, prediction_target_var])
    data$seasonally_differenced_log_prediction_target <- 
        ts(c(rep(NA, 52),
             data$log_prediction_target[seq(from = 53, to = nrow(data))] -
                 data$log_prediction_target[seq(from = 1, to = nrow(data) - 52)]),
           frequency = 52)
    
    
    ### make leave-one-season-out cross-validated predictions from SARIMA
    ## incidence bins for predictions of incidence in individual weeks and at the peak
    incidence_bins <- data.frame(
        lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
        upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf))
    
    ## at which weeks within the season do we form predictions?
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    
    ## data frame to describe predictions
    ## allocate more than enough space up front,
    ## delete extra later
    na_vec <- rep(NA,
                  length(unique(data$season)) *
                      (last_analysis_time_season_week - first_analysis_time_season_week + 1)
    )
    
    onset_bin_values <- c(40:52, 1:20, "none")
    onset_bin_names <- paste0("onset_week_", onset_bin_values)
    num_onset_bins <- length(onset_bin_names)
    
    peak_week_bin_values <- c(40:52, 1:20)
    peak_week_bin_names <- paste0("peak_week_", peak_week_bin_values)
    num_peak_week_bins <- length(peak_week_bin_names)
    
    incidence_bin_names <- paste0("incidence_", seq(from = 0, to = 13, by = 0.1))
    num_incidence_bins <- length(incidence_bin_names)
    
    predictions_df <- data.frame(
        model = "sarima",
        analysis_time_season = na_vec,
        analysis_time_season_week = na_vec,
        prediction_week_ph_1 = na_vec,
        prediction_week_ph_2 = na_vec,
        prediction_week_ph_3 = na_vec,
        prediction_week_ph_4 = na_vec,
        stringsAsFactors = FALSE) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_onset_bins)) %>%
                `colnames<-`(paste0("log_prob_", onset_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_peak_week_bins)) %>%
                `colnames<-`(paste0("log_prob_", peak_week_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_incidence_bins)) %>%
                `colnames<-`(paste0("log_prob_peak_inc_", incidence_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_incidence_bins)) %>%
                `colnames<-`(paste0("log_prob_ph1_", incidence_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_incidence_bins)) %>%
                `colnames<-`(paste0("log_prob_ph2_", incidence_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_incidence_bins)) %>%
                `colnames<-`(paste0("log_prob_ph3_", incidence_bin_names))
        ) %>%
        cbind(
            as.data.frame(matrix(NA, nrow = length(na_vec), ncol = num_incidence_bins)) %>%
                `colnames<-`(paste0("log_prob_ph4_", incidence_bin_names))
        )
    
    results_save_row <- 1L
    
    for(analysis_time_season in unique(data$season)) {
        message(paste("Region", reg_num, analysis_time_season, "::", Sys.time()))
        ## Only do something if there is something to predict in the season that would be held out
        if(!all(is.na(data$seasonally_differenced_log_prediction_target[data$season == analysis_time_season]))) {
            ## load SARIMA fit
            sarima_fit <- readRDS(file = paste0(
                filepath,
                "sarima-region",
                sprintf("%02d", reg_num),
                "-fit-leave-out-",
                gsub("/", "-", analysis_time_season),
                ".rds"))
            
            ## make predictions for each prediction target in the left-out season
            ## for each possible "last observed" week, starting with the last week of the previous season
            first_season_ind <- min(which(data$season == analysis_time_season))
            last_season_ind <- max(which(data$season == analysis_time_season))
            
            last_analysis_time_season_week_in_data <- max(data$season_week[data$season == analysis_time_season])
            for(analysis_time_season_week in seq(from = first_analysis_time_season_week, to = min(last_analysis_time_season_week, last_analysis_time_season_week_in_data) - 1)) {
                analysis_time_ind <- which(data$season == analysis_time_season &
                                               data$season_week == analysis_time_season_week)
                
                ## only do further work if the data ending at the analysis time contain
                ## enough recent observations to use the sarima fit.
                ## I'm not explicitly checking for seasonal lags...
                ## Subsetting a vector of class "ts" gives you a numeric instead of a "ts"
                ## Be sure to reset class to ts with frequency 52.
                new_target_data <- ts(data$seasonally_differenced_log_prediction_target[seq_len(analysis_time_ind)],
                                      frequency = 52)
                
                ar_order <- length(grep("^ar*", names(coef(sarima_fit))))
                
                #      !any(is.na(tail(new_target_data, ar_order)))
                ## keep track of if we made any predictions with this combination of season and week
                made_predictions <- FALSE
                
                if(!any(is.na(tail(new_target_data, ar_order)))) {
                    ## Update anova fit object with seasonally differenced data up
                    ## through analysis_time_ind
                    updated_sarima_fit <- Arima(
                        new_target_data,
                        model = sarima_fit)
                    
                    ## The update process via call to Arima can somehow result in 0/negative
                    ## variance estimates.  Reset to original values
                    updated_sarima_fit$sigma2 <- sarima_fit$sigma2
                    updated_sarima_fit$var.coef <- sarima_fit$var.coef
                    
                    ## simulate n_sims trajectories recursively from sarima
                    max_prediction_horizon <- max(4L,
                                                  last_analysis_time_season_week + 1 - analysis_time_season_week)
                    
                    raw_trajectory_samples <- sample_predictive_trajectories_arima(
                        updated_sarima_fit,
                        h = max_prediction_horizon,
                        npaths = n_sims)
                    
                    ## Sampled trajectories are of seasonally differenced log incidence
                    ## Get to trajectories for originally observed incidence ("inc") by
                    ## adding seasonal lag of incidence and exponentiating
                    inc_trajectory_samples <- raw_trajectory_samples
                    for(prediction_horizon in seq_len(max_prediction_horizon)) {
                        inc_trajectory_samples[, prediction_horizon] <- raw_trajectory_samples[, prediction_horizon] +
                            data$log_prediction_target[analysis_time_ind + prediction_horizon - 52]
                    }
                    inc_trajectory_samples <- exp(inc_trajectory_samples)
                    
                    
                    ## NA values could result if there were NAs in the data at seasonal lags that were used.
                    ## If that happened, forget about doing prediction for the given quantity
                    sample_inds_with_na <- apply(inc_trajectory_samples, 1, function(x) any(is.na(x)))
                    
                    ## Predictions for things about the whole season
                    if(!all(sample_inds_with_na)) {
                        made_predictions <- TRUE
                        
                        ## subset to sampled trajectories that are usable/do not have NAs
                        subset_inc_trajectory_samples <- inc_trajectory_samples[!sample_inds_with_na, ]
                        
                        ## Augment trajectory samples with previous observed incidence values
                        ## This is where we should be adding in some boostrapping to account
                        ## for backfill.
                        season_start_ind <- which(data$season == analysis_time_season &
                                                      data$season_week == 1)
                        if(season_start_ind < analysis_time_ind) {
                            subset_inc_trajectory_samples <- cbind(
                                matrix(
                                    rep(data[seq(from = season_start_ind, to = analysis_time_ind), prediction_target_var], each = n_sims),
                                    nrow = n_sims
                                ),
                                subset_inc_trajectory_samples
                            )
                        }
                        
                        ## Get onset week for each simulated trajectory
                        onset_week_by_sim_ind <- apply(subset_inc_trajectory_samples, 1, function(trajectory) {
                            get_onset_week(
                                incidence_trajectory = trajectory,
                                baseline = baselines$baseline[baselines$season=="2015/2016"], ## hard coded!
                                onset_length = 3L
                            )
                        })
                        ## Above gives onset week within season -- convert to week within year.
                        onset_week_by_sim_ind <- (data$week[data$season == analysis_time_season])[onset_week_by_sim_ind]
                        
                        ## Get peak week and height at peak week for each simulated trajectory
                        peak_week_by_sim_ind <- apply(subset_inc_trajectory_samples, 1, which.max)
                        # predictions_df[results_save_row, paste0("peak_week_", seq_len(n_sims))] <-
                        #   peak_week_by_sim_ind
                        
                        peak_week_height_by_sim_ind <- subset_inc_trajectory_samples[cbind(seq_len(n_sims), peak_week_by_sim_ind)]
                        # predictions_df[results_save_row, paste0("unbinned_peak_height_", seq_len(n_sims))] <-
                        #   peak_week_height_by_sim_ind
                        
                        peak_week_height_by_sim_ind <- sapply(peak_week_height_by_sim_ind,
                                                              function(height) {
                                                                  which(incidence_bins$lower <= height &
                                                                            incidence_bins$upper > height)
                                                              })
                        #          predictions_df[results_save_row, paste0("peak_height_", seq_len(n_sims))] <-
                        #            peak_week_height_by_sim_ind
                        
                        ## Above gives peak week within season -- convert to week within year.
                        peak_week_by_sim_ind <- (data$week[data$season == analysis_time_season])[peak_week_by_sim_ind]
                        
                        ## Get log of predictive probabilities for each bin
                        predictions_df[results_save_row, paste0("log_prob_", onset_bin_names)] <-
                            log(table(factor(onset_week_by_sim_ind, levels = onset_bin_values))) -
                            log(nrow(subset_inc_trajectory_samples))
                        predictions_df[results_save_row, paste0("log_prob_", peak_week_bin_names)] <-
                            log(table(factor(peak_week_by_sim_ind, levels = peak_week_bin_values))) -
                            log(nrow(subset_inc_trajectory_samples))
                        predictions_df[results_save_row, paste0("log_prob_peak_inc_", incidence_bin_names)] <-
                            log(table(factor(peak_week_height_by_sim_ind, levels = seq_len(num_incidence_bins)))) -
                            log(nrow(subset_inc_trajectory_samples))
                    }
                    
                    ## Predictions for incidence in an individual week at prediction horizon ph = 1, ..., 4
                    for(ph in 1:4) {
                        sample_inds_with_na <- is.na(inc_trajectory_samples[, ph])
                        
                        if(!all(sample_inds_with_na)) {
                            made_predictions <- TRUE
                            
                            ## get sampled incidence values at prediction horizon that are usable/not NAs
                            ph_inc_samples <- inc_trajectory_samples[!sample_inds_with_na, ph]
                            ph_inc_bin_samples <- sapply(ph_inc_samples, function(inc) {
                                which(
                                    incidence_bins$lower <= inc &
                                        incidence_bins$upper > inc
                                )
                            })
                            
                            ## Get log of predictive probabilities for each bin
                            predictions_df[results_save_row, paste0("log_prob_ph", ph, "_", incidence_bin_names)] <-
                                log(table(factor(ph_inc_bin_samples, levels = seq_len(num_incidence_bins)))) -
                                log(sum(!sample_inds_with_na))
                        }
                    } # ph loop
                }
                
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
        }
    } # analysis_time_season
    
    ## if there are extra rows in the predictions_df, delete them
    if(results_save_row <= nrow(predictions_df)) {
        predictions_df <- predictions_df[
            -seq(from = results_save_row, to = nrow(predictions_df)),
            ,
            drop = FALSE
            ]
    }
    
    sarima_predictions_df <- predictions_df
    filename <- paste0(
        "inst/estimation/loso-predictions/sarima-region",
        sprintf("%02d", reg_num),
        "-loso-predictions.rds")
    saveRDS(sarima_predictions_df, file = filename)
    
}
