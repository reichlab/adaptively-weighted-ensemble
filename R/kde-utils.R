## code for KDE prediction models
## Nicholas Reich
## 7 October 2016


## General strategy

## ESTIMATION
## write a function 
##   inputs: a dataset
##   outputs: a list of KDE fits, one for each target
## iterate function across all regions and years to create a set of LOSO fits
## 
## notes by target:
##     weekly incidence: use observations from target week +/- 1 week 
##     peak incidence: truncate above previous observations?
##     onset and peak week: discrete, but use continuous KDEs and round
##     onset week: there is a none category, consider including a point-mass based on past seasons not above baseline
##     
## PREDICTION
## write a predict function 
##   input: a fitted object from above 
##   input: a time at which to make predictions for (data not needed!)
##   input: n_sims
##   output: data_frame in format as others, specifying predictive distribution
## iterate across regions, left-out years to create predictions
## 


##' Fits and saves KDE models for all region-year combinations
##' 
##' @param data cleaned usflu dataset
##' @param region character string identifying a region
##' @param first_test_year along with first_test_week, defines the first time at which data is left out
##' @param first_test_week 
##' @param path filepath for saving 
##' 
##' @return nothing, just saving files
##' 
##' @export
##' 
fit_region_kdes <- function(data, region, first_test_year, first_test_week, path) {
    require(MMWRweek)
    require(dplyr)
    
    ### get proportion of region-seasons (across all regions and all seasons
    ### before first test year/week) with no onset
    onsets_by_region_season <- NULL
    for(region_val in unique(data$region)) {
        data_subset <- data[which(data$region == region_val),]
        ## subset to include only times before first test season and week
        idx <- which(data_subset$year == first_test_year & data_subset$week == first_test_week)
        data_subset <- data_subset[seq_len(idx - 1), , drop = FALSE]
        
        observed_seasonal_quantities_by_season <- t(sapply(
            as.character(unique(data_subset$season)),
            function(season) {
                get_observed_seasonal_quantities(
                    data = data_subset,
                    season = season,
                    first_CDC_season_week = 10,
                    last_CDC_season_week = 41,
                    onset_baseline = 
                        get_onset_baseline(region = region_val, season = season),
                    incidence_var = "weighted_ili",
                    incidence_bins = data.frame(
                        lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
                        upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
                    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
                )
            }
        ))
        
        onsets_by_region_season <- bind_rows(onsets_by_region_season,
            data.frame(
                region = region_val,
                season = rownames(observed_seasonal_quantities_by_season),
                onset = as.character(unname(unlist(observed_seasonal_quantities_by_season[, 1]))),
                stringsAsFactors = FALSE
            )
        )
    }
    
    ### subsetting data
    ## subset to region of interest
    dat <- data[which(data$region == region),]
    ## subset to include only times before first test season and week
    idx <- which(dat$year == first_test_year & dat$week == first_test_week)
    dat <- dat[seq_len(idx - 1), , drop = FALSE]

    ## assumes region is either "X" or "Region k" format
    reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
    
    ### loop over all seasons, including first season of test phase
    ### (i.e., create a fit that doesn't leave any training data out)
    for(season_left_out in c(unique(as.character(dat$season)), paste0(first_test_year, "/", first_test_year + 1))) {
        tmpdat <- dat
        
        ### create filename for saving and check to see if fit exists already
        filename <- paste0(
            path,
            "kde-",
            reg_string,
            "-fit-leave-out-",
            gsub("/", "-", season_left_out),
            ".rds")
        if(file.exists(filename)){
            message(paste(region, "leaving out", season_left_out, "already fit ::", Sys.time()))
            next
        }
        
        ### drop left-out season
        idx_to_drop <- tmpdat$season == season_left_out
        if(sum(idx_to_drop) > 0) {
          tmpdat <- tmpdat[-idx_to_drop,]
        }
        
        ### create fits
        kde_onset_week <- fit_kde_onset_week(tmpdat,
          prob_no_onset = mean(onsets_by_region_season$onset[onsets_by_region_season$season != season_left_out] == "none"))
        kde_peak_week <- fit_kde_peak_week(tmpdat)
        kde_log_peak_week_inc <- fit_kde_log_peak_week_inc(tmpdat)
        kde_log_weekly_inc <- fit_kde_weekly_inc(tmpdat)
            
        kde_fits <- list(onset_week = kde_onset_week, 
                         peak_week = kde_peak_week, 
                         log_peak_week_inc = kde_log_peak_week_inc, 
                         log_weekly_inc = kde_log_weekly_inc)
        
        ### save fits
        saveRDS(kde_fits, file = filename)
    }
}


#' Estimate KDE for onset week
#'
#' @param data data with one season left out
#' @param prob_no_onset probability to assign to "no onset" bin
#'
#' @return a fit from density(), along with points to evaluate at
fit_kde_onset_week <- function(data, prob_no_onset) {
    require(dplyr)
    ### get onset
    observed_seasonal_quantities_by_season <- t(sapply(
        unique(data$season),
        function(season) {
            get_observed_seasonal_quantities(
                data = data,
                season = season,
                first_CDC_season_week = 10,
                last_CDC_season_week = 41,
                onset_baseline = 
                    get_onset_baseline(region = data$region[1], season = season),
                incidence_var = "weighted_ili",
                incidence_bins = data.frame(
                    lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
                    upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
                incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
            )
        }
    ))
    
    ## vector of onset weeks
    onset_week <- as.numeric(unlist(observed_seasonal_quantities_by_season[, 1]))
    
    ### calculate density based on non-NA values
    onset_week <- onset_week[!is.na(onset_week)]
    return(list(kde=density(onset_week, na.rm=TRUE), prob_no_onset=prob_no_onset, x=onset_week))
}


#' Estimate KDE for peak week
#'
#' @param data 
#'
#' @return a fit from density(), along with points to evaluate at 
fit_kde_peak_week <- function(data) {
    require(dplyr)
    ### get peaks
    peak_week <- data %>% group_by(season) %>%
        ## filter to ensure estimated peaks are within CDC-specified forecasting season
        filter(season_week >= 10 & season_week <= 42) %>%
        summarize(peak_week = season_week[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        .$peak_week

    ### calculate density
    return(list(kde=density(peak_week, na.rm=TRUE), x=peak_week))
}

#' Estimate KDE for LOG peak week incidence
#'
#' @param data 
#'
#' @return a fit from density(), along with points to evaluate at
fit_kde_log_peak_week_inc <- function(data) {
    require(dplyr)
    ### get peaks
    peak_week_inc <- data %>% group_by(season) %>%
        summarize(peak_week_inc = weighted_ili[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        .$peak_week_inc
    
    ### calculate density
    return(list(kde=density(log(peak_week_inc), na.rm=TRUE), x=log(peak_week_inc)))
}

#' Fit a GAM for weekly incidence
#'
#' @param data 
#'
#' @return a mgcv smooth spline fit
fit_kde_weekly_inc <- function(data) {
    require(mgcv)
    
    ### fit model
    fm <- gam(log(weighted_ili) ~ s(season_week, bs = "cc"),
              data=data)
    
    return(fm)
}



#' Wrapper for prediction and evaluation of KDE fits
#'
#' @param data data, to be used for evaluation
#' @param region region of focus
#' @param path filepath to fitted models
#' @param n_sim number of simulations to run for predictive distributions
#'
#' @return NULL just saves a file
#'
#' @description The function follows this outline
#' \itemize{
#'       \item{determine set of years for which LOSO fits are available}
#'       \item{For those years, generate predictions made at each season-week. 
#'             These will be the same predictions repeated, as models not dependent on any inputs.}
#'       \item{Compare predictions to the real data}
#'       \item{Save object with comparisons}
#' }
#' 
#' @export
predict_region_kde <- function(data, region, path, n_sim) {
    
    require(dplyr)
    
    ### SETUP
    ## load baselines for onset calculations
    baselines <- read.csv(file="data-raw/cdc-baselines.csv")
    ## NOTE: assumes region is either "X" or "Region k" format
    reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
    idx <- which(baselines$region==reg_string & baselines$season=="2015/2016")
    reg_baseline <- baselines[idx, "baseline"]
    
    ## subset data to be only the region of interest
    dat <- data[which(data$region == region),]
    
    ### SET GLOBAL VALUES

    ## incidence bins for predictions of incidence in individual weeks and at the peak
    inc_bins <- c(0, seq(from = .05, to = 12.95, by = 0.1), Inf)
    incidence_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    ## for 2016-2017 season, first predictions due on 11/7/2016 (EW45 == SW15)
    ## using data posted on 11/4/2016 that includes up to EW43 == SW13
    ## last predictions due on 5/15/2017 (EW20 == SW 42)
    ## last predictions use data through EW40 == SW18
    ## first_analysis_time_season_week could be set to 13, but padding at front end
    first_analysis_time_season_week <- 10 # == week 40 of year
    last_analysis_time_season_week <- 41 # == week 19 or 18, depending on whether a 53 week season
    
    ## subset of season weeks for CDC competition
    cdc_season_weeks <- as.character(10:42)
    
    ### determine which years we have fits for
    fnames <- system(paste0('ls ', path), intern=TRUE)
    fnames_this_reg <- fnames[grep(paste0(reg_string,"-"), fnames)] ## need hyphen for Reg 1 vs. 10
    where_is_the_dot <- regexpr(pattern = "\\.", fnames_this_reg)
    seasons_this_reg <- substr(fnames_this_reg, start=where_is_the_dot - 9, stop=where_is_the_dot-1)
    seasons_this_reg <- gsub("-", "/", seasons_this_reg)
    
    ### calculate observed values of targets across seasons for later comparison
    ## season week is week of flu season (week 1 is in Oct)
    ## week is epidemic week of calendar year (week 1 is in January)
    
    ## onset week
    onset_weeks <- dat %>% group_by(season) %>% 
        mutate(ili_lag1 = lag(weighted_ili, 1),
               ili_lag2 = lag(weighted_ili, 2),
               ili_lag3 = lag(weighted_ili, 3),
               onset = ili_lag1 > reg_baseline & ili_lag2 > reg_baseline & ili_lag3 > reg_baseline) %>%
        filter(onset) %>%
        summarise(season_week = min(season_week)-2) %>% ## may not do the right thing if onset is in first 2 weeks
        left_join(dat) %>%
        select(season, season_week, week)
    
    ## peak week
    peak_weeks <- dat %>% group_by(season) %>%
        summarize(season_week = season_week[which(weighted_ili == max(weighted_ili, na.rm=TRUE))]) %>%
        left_join(dat) %>%
        select(season, season_week, week)
    
    ## peak week incidence
    peak_week_inc <- dat %>% group_by(season) %>%
        summarize(peak_week_inc = weighted_ili[which(weighted_ili == max(weighted_ili, na.rm=TRUE))])
    
    ### generate predictions made at each season-week
    ## make dataset for storage
    predictions_df <- make_predictions_dataframe(dat, 
                                                 model_name="kde", 
                                                 incidence_bin_names=incidence_bin_names)
    
    results_save_row <- 1L
    
    for(analysis_time_season in seasons_this_reg) {
        message(paste(region, analysis_time_season, "::", Sys.time()))
        ## Only do something if there is something to predict in the season that would be held out
        if(!all(is.na(dat$weighted_ili[dat$season == analysis_time_season]))) {
            ## load KDE fit
            kde_fit <- readRDS(file = paste0(
                path,
                "kde-",
                reg_string,
                "-fit-leave-out-",
                gsub("/", "-", analysis_time_season),
                ".rds"))
            
            ## make predictions for each prediction target in the left-out season
            ## for each possible "last observed" week, starting with the last week of the previous season
            first_season_ind <- min(which(dat$season == analysis_time_season))
            last_season_ind <- max(which(dat$season == analysis_time_season))
            
            last_analysis_time_season_week_in_dat <- max(dat$season_week[dat$season == analysis_time_season])
            
            ## make prediction for this analysis_time_season
            onset_week_preds <- predict_kde_onset_week(kde_fit$onset_week, n_sim)
            peak_week_preds <- predict_kde_peak_week(kde_fit$peak_week, n_sim)
            peak_week_inc_preds <- predict_kde_log_peak_week_inc(kde_fit$log_peak_week_inc, 
                                                                 bins=inc_bins,
                                                                 bin_names=incidence_bin_names, 
                                                                 n_sim)
            weekly_inc_preds <- predict_kde_log_weekly_inc(fm = kde_fit$log_weekly_inc, 
                                                           season_weeks = 1:53, 
                                                           bins = inc_bins, 
                                                           bin_names = incidence_bin_names,
                                                           n_sim = n_sim)
            
            ## calculate observed values to compare against
            observed_onset_week <- unlist(onset_weeks[which(onset_weeks$season==analysis_time_season),"season_week"])
            observed_peak_week <- unlist(peak_weeks[which(peak_weeks$season==analysis_time_season),"season_week"])
            observed_peak_week_inc <- unlist(peak_week_inc[which(peak_week_inc$season==analysis_time_season),"peak_week_inc"])
            
            ### calculate log score for unchanging predictions
            
            ## calculate onset week log scores
            if(length(observed_onset_week)==0) { ## if no observed onset, assign none
                obs_onset_week_char <- "none"
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            } else if(observed_onset_week<10) { ## if onset prior to week 10, assign log-score of NA
                log_score_for_onset_week <- NA
            } else if(observed_onset_week>42) { ## if onset after week 42, assign "none"
                obs_onset_week_char <- "none"
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            } else { ## otherwise, onset week in [10, 42]
                obs_onset_week_char <- as.character(observed_onset_week)
                log_score_for_onset_week <- compute_competition_log_score(log(onset_week_preds[c(cdc_season_weeks, "none")]),
                                                                          observed_bin = obs_onset_week_char,
                                                                          prediction_target = "onset_week")
            }

            ## calculate peak week log scores
            ## assign the index based on output from predict_kde_peak 
            log_score_for_peak_week <- ifelse(observed_peak_week<10 | observed_peak_week > 42,
                                              NA,
                                              compute_competition_log_score(log(peak_week_preds[cdc_season_weeks]),
                                                                     observed_bin = as.character(observed_peak_week),
                                                                     prediction_target = "peak_week"))
            
            ## calculate peak week incidence log scores
            peak_inc_bin_char <- get_inc_bin(observed_peak_week_inc)
            log_score_for_peak_week_inc <- compute_competition_log_score(log(peak_week_inc_preds),
                                                                         observed_bin = peak_inc_bin_char,
                                                                         prediction_target="peak_inc")

            ## sequence to loop over
            time_seq <- seq(from = first_analysis_time_season_week, to = min(last_analysis_time_season_week, last_analysis_time_season_week_in_dat) - 1)
            
            for(analysis_time_season_week in time_seq) {
                ## assign log scores
                predictions_df[results_save_row, "onset_log_score"] <- log_score_for_onset_week
                predictions_df[results_save_row, "peak_week_log_score"] <- log_score_for_peak_week
                predictions_df[results_save_row, "peak_inc_log_score"] <- log_score_for_peak_week_inc
                
                ## calculate and assign peak week incidence log scores
                analysis_time_ind <- which(dat$season == analysis_time_season &
                                               dat$season_week == analysis_time_season_week)

                ## 1 week ahead
                observed_ph_1_inc <- dat[analysis_time_ind+1, "weighted_ili"]
                ph_1_inc_bin <- get_inc_bin(observed_ph_1_inc)
                ph_1_season_week <- analysis_time_season_week+1
                log_score_for_ph_1_inc <- ifelse(is.na(ph_1_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_1_season_week]),
                                                                        observed_bin = ph_1_inc_bin,
                                                                        prediction_target="ph1_inc"))
                predictions_df[results_save_row, "ph_1_inc_log_score"] <- log_score_for_ph_1_inc
                
                ## 2 week ahead
                observed_ph_2_inc <- dat[analysis_time_ind+2, "weighted_ili"]
                ph_2_inc_bin <- get_inc_bin(observed_ph_2_inc)
                ph_2_season_week <- analysis_time_season_week+2
                log_score_for_ph_2_inc <- ifelse(is.na(ph_2_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_2_season_week]),
                                                                        observed_bin = ph_2_inc_bin,
                                                                        prediction_target="ph2_inc"))
                predictions_df[results_save_row, "ph_2_inc_log_score"] <- log_score_for_ph_2_inc
                
                ## 3 week ahead
                observed_ph_3_inc <- dat[analysis_time_ind+3, "weighted_ili"]
                ph_3_inc_bin <- get_inc_bin(observed_ph_3_inc)
                ph_3_season_week <- analysis_time_season_week+3
                log_score_for_ph_3_inc <- ifelse(is.na(ph_3_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_3_season_week]),
                                                                        observed_bin = ph_3_inc_bin,
                                                                        prediction_target="ph3_inc"))
                predictions_df[results_save_row, "ph_3_inc_log_score"] <- log_score_for_ph_3_inc
                
                ## 4 week ahead
                observed_ph_4_inc <- dat[analysis_time_ind+4, "weighted_ili"]
                ph_4_inc_bin <- get_inc_bin(observed_ph_4_inc)
                ph_4_season_week <- analysis_time_season_week+4
                log_score_for_ph_4_inc <- ifelse(is.na(ph_4_inc_bin), ## check NA
                                                 NA,
                                                 compute_competition_log_score(log(weekly_inc_preds[, ph_4_season_week]),
                                                                        observed_bin = ph_4_inc_bin,
                                                                        prediction_target="ph4_inc"))
                predictions_df[results_save_row, "ph_4_inc_log_score"] <- log_score_for_ph_4_inc
                
                ## set other fixed stuff
                predictions_df[results_save_row, "analysis_time_season"] <- analysis_time_season
                predictions_df[results_save_row, "analysis_time_season_week"] <- analysis_time_season_week
                predictions_df[results_save_row, "prediction_week_ph_1"] <- analysis_time_season_week + 1
                predictions_df[results_save_row, "prediction_week_ph_2"] <- analysis_time_season_week + 2
                predictions_df[results_save_row, "prediction_week_ph_3"] <- analysis_time_season_week + 3
                predictions_df[results_save_row, "prediction_week_ph_4"] <- analysis_time_season_week + 4
                
                results_save_row <- results_save_row + 1
                
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
    
    ### save and return
    kde_predictions_df <- predictions_df
    filename <- paste0(
        "inst/estimation/loso-predictions/kde-",
        reg_string,
        "-loso-predictions.rds")
    saveRDS(kde_predictions_df, file = filename)
    
}

#' Get log scores and full predictive distributions for each prediction target
#' 
#' for a given season using a predictive method that works by directly 
#' simulating predictive distributions of each target.  Results are
#' stored in a data frame, saved in a .rds file with a name like
#' "model_name-region-season-loso-predictions.rds"
#' Results have columns indicating the analysis time season and season week,
#' model name, log scores for each prediction target, the "log score" used
#' in the competition (adding probabilities from adjacent bins) for each
#' prediction target, as well as the log of the probability assigned to each
#' bin.
#' 
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
#'   "Region k" where k in {1, ..., 10}
#' @param prediction_target_var string specifying the name of the variable in
#'   data for which we want to make predictions
#' @param incidence_bins a data frame with variables lower and upper defining
#'   lower and upper endpoints to use in binning incidence
#' @param incidence_bin_names a character vector with a name for each incidence
#'   bin
#' @param n_sims integer number of samples to simulate
#' @param model_name name of model, stored in the results data frame
#' @param fits_path path to directory where fitted models are stored
#' @param prediction_save_path path to directory where results will be saved
#' 
#' @return none
#' 
#' @export
get_log_scores_via_direct_simulation <- function(
    analysis_time_season,
    first_analysis_time_season_week = 10, # == week 40 of year
    last_analysis_time_season_week = 41, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
    region,
    prediction_target_var,
    incidence_bins,
    incidence_bin_names,
    n_sims,
    model_name,
    fits_path,
    prediction_save_path
) {
    require(dplyr)
    
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
    
    ## Only do something if there is something to predict in the season that would be held out
    if(!all(is.na(data[data$season == analysis_time_season, prediction_target_var]))) {
        ## load KDE fit
        ## NOTE: assumes region is either "X" or "Region k" format
        reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
        kde_fit <- readRDS(file = paste0(
            fits_path,
            "kde-",
            reg_string,
            "-fit-leave-out-",
            "2011-2012", # always use 2011/2012 season fit -- based on all training data
#            gsub("/", "-", analysis_time_season),
            ".rds"))

        ### calculate log score for predictions that don't change depending on time of year
        
        ## Get predictions and log scores for onset
        onset_week_bins <- c(as.character(10:42), "none")
        onset_week_preds <- predict_kde_onset_week(kde_fit$onset_week, n_sims) # returns probs for weeks 1:52
        onset_bin_log_probs <- log(onset_week_preds[onset_week_bins]) # subset to the weeks we care about (should maybe sum probs for weeks before 10?)
        onset_bin_log_probs <- onset_bin_log_probs - logspace_sum(onset_bin_log_probs) # re-normalize after subsetting
        predictions_df[, paste0("onset_bin_", onset_week_bins, "_log_prob")] <-
            rep(onset_bin_log_probs, each = nrow(predictions_df))
        predictions_df[, "onset_log_score"] <-
            onset_bin_log_probs[ as.character(observed_seasonal_quantities$observed_onset_week) ]
        predictions_df[, "onset_competition_log_score"] <-
            compute_competition_log_score(onset_bin_log_probs,
                                          as.character(observed_seasonal_quantities$observed_onset_week),
                                          "onset_week")
       
        ## Get log scores for peak week
        peak_week_bins <- as.character(10:42)
        peak_week_preds <- predict_kde_peak_week(kde_fit$peak_week, n_sims) ## returns probs for all weeks
        peak_week_bin_log_probs <- log(peak_week_preds[peak_week_bins]) ## subset to only competition weeks
        peak_week_bin_log_probs <- peak_week_bin_log_probs - logspace_sum(peak_week_bin_log_probs) # re-normalize after subsetting
        names(peak_week_bin_log_probs) <- as.character(peak_week_bins)
        predictions_df[, paste0("peak_week_bin_", peak_week_bins, "_log_prob")] <-
            rep(peak_week_bin_log_probs, each = nrow(predictions_df))
        predictions_df[, "peak_week_log_score"] <-
            logspace_sum(peak_week_bin_log_probs[ as.character(observed_seasonal_quantities$observed_peak_week)] )
        predictions_df[, "peak_week_competition_log_score"] <-
            compute_competition_log_score(peak_week_bin_log_probs,
                                          as.character(observed_seasonal_quantities$observed_peak_week),
                                          "peak_week")
        
        ## Get log scores for peak week incidence
        peak_week_inc_preds <- predict_kde_log_peak_week_inc(kde_fit$log_peak_week_inc, 
                                                             bins=c(0, incidence_bins$upper),
                                                             bin_names=incidence_bin_names, 
                                                             n_sims)
        peak_inc_bin_log_probs <- log(peak_week_inc_preds)
        predictions_df[, paste0("peak_inc_bin_", incidence_bin_names, "_log_prob")] <-
            rep(peak_inc_bin_log_probs, each = nrow(predictions_df))
        predictions_df[, "peak_inc_log_score"] <-
            peak_inc_bin_log_probs[
                as.character(observed_seasonal_quantities$observed_peak_inc_bin)]
        predictions_df[, "peak_inc_competition_log_score"] <-
            compute_competition_log_score(peak_inc_bin_log_probs,
                                          observed_seasonal_quantities$observed_peak_inc_bin,
                                          "peak_inc")
    
        ## make prediction for this analysis_time_season
        weekly_inc_preds <- predict_kde_log_weekly_inc(fm = kde_fit$log_weekly_inc, 
                                                       season_weeks = 1:53, 
                                                       bins = c(0, incidence_bins$upper), 
                                                       bin_names = incidence_bin_names,
                                                       n_sim = n_sims)
        
        
        ## figure out weeks for looping
        last_analysis_time_season_week_in_data <- max(data$season_week[data$season == analysis_time_season])
        time_seq <- seq(from = first_analysis_time_season_week, 
                        to = min(last_analysis_time_season_week, last_analysis_time_season_week_in_data) - 1)
        
        ## loop over all times in season
        results_save_row <- 1L
        for(analysis_time_season_week in time_seq) {

            ## calculate current time 
            analysis_time_ind <- which(data$season == analysis_time_season &
                                           data$season_week == analysis_time_season_week)
            
            ## Predictions for incidence in an individual week at prediction horizon ph = 1, ..., 4
            for(ph in 1:4) {
                ## get observed value/bin
                observed_ph_inc <- data[analysis_time_ind + ph, prediction_target_var]
                observed_ph_inc_bin <- get_inc_bin(observed_ph_inc, return_character = TRUE)
                
                if(!is.na(observed_ph_inc)) {

                    ## get sampled incidence values at prediction horizon that are usable/not NAs
                    pred_week_idx <- analysis_time_season_week + ph
                    ph_inc_bin_preds <- weekly_inc_preds[, pred_week_idx]
                    ## removed call to get_inc_bin() b/c predict_kde_log_weekly_inc() does it for us.

                    ## get log score
                    ph_inc_bin_log_probs <- log(ph_inc_bin_preds)
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
         
            ## set other fixed stuff
            predictions_df[results_save_row, "analysis_time_season"] <- analysis_time_season
            predictions_df[results_save_row, "analysis_time_season_week"] <- analysis_time_season_week
            predictions_df[results_save_row, "prediction_week_ph_1"] <- analysis_time_season_week + 1
            predictions_df[results_save_row, "prediction_week_ph_2"] <- analysis_time_season_week + 2
            predictions_df[results_save_row, "prediction_week_ph_3"] <- analysis_time_season_week + 3
            predictions_df[results_save_row, "prediction_week_ph_4"] <- analysis_time_season_week + 4
            
            results_save_row <- results_save_row + 1
            
        } # analysis_time_season_week
    }
   
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
                          model_name, "-", region_str, "-", season_str, "-loso-predictions.rds"))
}

