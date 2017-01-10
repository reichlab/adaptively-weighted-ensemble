## function to make current predictions from KDE
## Nicholas Reich
## started: November 2016

#' Wrapper for prediction of a single time-point from KDE fits
#'
#' @param path filepath to fitted models
#' @param season character string of season for current fit, in format "XXXX/YYYY"
#' @param season_week numeric season week of first week for which predictions are desired
#' @param n_sim number of simulations to run for predictive distributions
#'
#' @return NULL just saves a file
#'
make_one_kde_prediction_file <- function(path, season, season_week, n_sim){
    ## set globals for incidence bins and names
    inc_bins <- c(0, seq(from = .05, to = 12.95, by = 0.1), Inf)
    incidence_bin_names <- as.character(seq(from = 0, to = 13, by = 0.1))
    
    region_strings <- c("National", paste0("Region", 1:10))
    cdc_region_strings <- c("US National", paste("HHS Region", 1:10))
    
    epiweek <- (season_week-1+30)%%52 + 1
    fname <- paste0(path,"/EW", epiweek, "-KOTkde-", Sys.Date(), ".csv")
    
    for(i in 1:length(region_strings)) {
        region <- region_strings[i]
        cdc_region <- cdc_region_strings[i]
        ## load KDE fit
        kde_fit <- readRDS(file = paste0(
            path,
            "kde-",
            region,
            "-fit-leave-out-",
            gsub("/", "-", season),
            ".rds"))
        tmp <- make_region_prediction(kde_fit, 
            cdc_region_string=cdc_region, 
            inc_bins=inc_bins, 
            inc_bin_names=incidence_bin_names, 
            season_week = season_week,
            n_sim=n_sim)
        add_colnames <- ifelse(i==1, TRUE, FALSE) ## only include rownames in first file
        write.table(tmp, file=fname, quote=FALSE, sep=",", row.names = FALSE,
            append=!add_colnames, col.names=add_colnames)
    }
}

make_region_prediction <- function(fit, cdc_region_string, inc_bins, inc_bin_names, season_week, n_sim) {
    ## read in output template
    out <- read.csv("data-raw/region-prediction-template.csv")
    out$Location <- cdc_region_string
    
    ### ONSETS
    onset_week_preds <- predict_kde_onset_week(fit$onset_week, n_sim)
    ## truncating to weeks observed and adding delta to ensure no -Inf log scores
    onset_week_preds_used <- onset_week_preds[c(10:42, 53)] + 1/n_sim
    onset_week_preds_std <- onset_week_preds_used/sum(onset_week_preds_used)
    idx_onset_bins <- which(out$Target=="Season onset" & out$Type=="Bin")
    out[idx_onset_bins, "Value"] <- onset_week_preds_std
    idx_onset_point <- which(out$Target=="Season onset" & out$Type=="Point")
    onset_season_week <- calc_median_from_binned_probs(onset_week_preds_std)
    out[idx_onset_point, "Value"] <- (onset_season_week-1+30)%%52 +1 ## dependent on 52 week season
    
    ### PEAKS
    peak_week_preds <- predict_kde_peak_week(fit$peak_week, n_sim)
    ## truncating to weeks observed and adding delta to ensure no -Inf log scores
    peak_week_preds_used <- peak_week_preds[10:42] + 1/n_sim
    peak_week_preds_std <- peak_week_preds_used/sum(peak_week_preds_used)
    idx_peak_bins <- which(out$Target=="Season peak week" & out$Type=="Bin")
    out[idx_peak_bins, "Value"] <- peak_week_preds_std
    idx_peak_point <- which(out$Target=="Season peak week" & out$Type=="Point")
    peak_season_week <- calc_median_from_binned_probs(peak_week_preds_std)
    out[idx_peak_point, "Value"] <- (peak_season_week-1+30)%%52 +1 ## dependent on 52 week season
    
    ### PEAK INC
    peak_inc_preds <- predict_kde_log_peak_week_inc(fit$log_peak_week_inc, 
        bins=inc_bins,
        bin_names=inc_bin_names, 
        n_sim)
    peak_inc_preds_used <- peak_inc_preds + 1/n_sim
    peak_inc_preds_std <- peak_inc_preds_used/sum(peak_inc_preds_used)
    idx_peak_inc_bins <- which(out$Target=="Season peak percentage" & out$Type=="Bin")
    ## check that names match
    if(any(names(peak_inc_preds) != out[idx_peak_inc_bins, "Bin_start_incl"])){
        error("Peak incidencebin names don't match.")
    }
    out[idx_peak_inc_bins, "Value"] <- peak_inc_preds_std
    idx_peak_inc_point <- which(out$Target=="Season peak percentage" & out$Type=="Point")
    out[idx_peak_inc_point, "Value"] <- calc_median_from_binned_probs(peak_inc_preds_std)
    
    
    ## WEEKLY INCIDENCE
    weekly_inc_preds <- predict_kde_log_weekly_inc(fm = fit$log_weekly_inc, 
                                                   season_weeks = season_week:(season_week+3), 
                                                   bins = inc_bins, 
                                                   bin_names = inc_bin_names,
                                                   n_sim = n_sim)
    weekly_inc_preds <- weekly_inc_preds + 1/n_sim
    weekly_inc_preds <- apply(weekly_inc_preds, 
        FUN=function(x) x/sum(x),
        MAR=2)
    for(i in 1:4){
        idx_weekly_bins <- which(out$Target==paste(i, "wk ahead") & out$Type=="Bin")
        out[idx_weekly_bins, "Value"] <- weekly_inc_preds[,i]
        idx_weekly_inc_point <- which(out$Target==paste(i, "wk ahead") & out$Type=="Point")
        weekly_inc_preds_vec <- weekly_inc_preds[,i]
        names(weekly_inc_preds_vec) <- inc_bin_names
        out[idx_weekly_inc_point, "Value"] <- calc_median_from_binned_probs(weekly_inc_preds_vec)
    }
    
    return(out)
}

