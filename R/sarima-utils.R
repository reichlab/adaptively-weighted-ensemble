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
    
    if(as.integer(substr(analysis_time_season, 1, 4)) >=
        as.integer(substr(params$first_test_season, 1, 4))) {
        ## if we're doing prediction during the test period,
        ## get fit obtained using all training data
        fit_season_name <- "none"
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
