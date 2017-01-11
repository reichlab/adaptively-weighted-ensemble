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
    
    ### subsetting data
    ## subset to region of interest
    dat <- data[which(data$region == region),]
    ## subset to include only times before first test season and week
    idx <- which(as.Date(as.character(dat$time)) <= MMWRweek2Date(first_test_year, first_test_week))
    dat <- dat[idx,]
    
    ## assumes region is either "X" or "Region k" format
    reg_string <- ifelse(region=="X", "National", gsub(" ", "", region))
    
    ### loop over all seasons
    for(season_left_out in unique(dat$season)) {
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
        tmpdat <- tmpdat[-idx_to_drop,]
        
        ### create fits
        kde_onset_week <- fit_kde_onset_week(tmpdat)
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
#'
#' @return a fit from density(), along with points to evaluate at
fit_kde_onset_week <- function(data) {
    require(dplyr)
    ### get onset
    onsets <- data %>% group_by(season) %>%
        mutate(baseline = get_onset_baseline(region = region, season = season),
               ili_lag1 = lag(weighted_ili, 1),
               ili_lag2 = lag(weighted_ili, 2),
               ili_lag3 = lag(weighted_ili, 3)) %>%
        ## filter to truncate onsets to only occur between 10 and 42
        filter(season_week <= 44 & season_week >= 12) %>% 
        mutate(onset = ili_lag1 > baseline & ili_lag2 > baseline & ili_lag3 > baseline) 
    
    ## vector of onset weeks
    onset_week <- onsets %>%
        filter(onset) %>% 
        summarize(first_week = min(season_week)-2) %>% ## may not do the right thing if onset is in first 2 weeks
        .$first_week
    
    ## also, calculate the probability of no onset
    any_onset <- onsets %>% 
        summarize(onset_occurred = any(onset)) %>%
        ungroup() 
    prob_no_onset <- 1-sum(any_onset$onset_occurred, na.rm=TRUE)/nrow(any_onset)

    ### calculate density
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


###########################
## BEGIN PREDICTION CODE ##
###########################

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
