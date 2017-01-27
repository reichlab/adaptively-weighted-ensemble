#' Functions for predicting from KDE fits
#'


#' Compute predictive distribution for onset_week from KDE fit
#'
#' @param fm output from fit_kde_onset_week()
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of 53 probabilities representing the predictive distribution of onset week, 
#'         in bins representing season week 1 through 52 (or EW31 through EW30 -- or EW29 if a 53 week year) 
#'         with the final entry representing the probability of no onset
#'
predict_kde_onset_week <- function(fm, n_sim) {
    kde_onset_week <- fm$kde
    prob_no_onset <- fm$prob_no_onset
    onsets <- fm$x
    x.new <- rnorm(n_sim, 
                   mean=sample(onsets, size = n_sim, replace = TRUE), 
                   sd=kde_onset_week$bw)
    pred_onsets_rounded <- round(x.new)

    ## removing early/late onset draws
    idx_out_of_bounds <- which(pred_onsets_rounded<=0 | pred_onsets_rounded>52)
    if(length(idx_out_of_bounds)>0){
      pred_onsets_rounded <- pred_onsets_rounded[-idx_out_of_bounds] 
    }
    
    onsets_binned <- tabulate(pred_onsets_rounded, nbins=52)
    ## calculate number of "nones" to pad for correct probability
    nones_to_add <- round(prob_no_onset/(1-prob_no_onset)*length(pred_onsets_rounded))

    onset_bin_probs <- c(onsets_binned, nones_to_add)/(nones_to_add+length(pred_onsets_rounded))
    names(onset_bin_probs) <- c(1:52, "none")
    
    return(onset_bin_probs)    
}

#' Compute predictive distribution for peak week from KDE fit
#'
#' @param fm output from fit_kde_peak_week()
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of 52 probabilities representing the predictive distribution of peak week, 
#'         in bins representing season week 1 through 52 (or EW31 through EW30 -- or EW29 if a 53 week year) 
#'         with the final entry representing the probability of no onset#'
predict_kde_peak_week <- function(fm, n_sim) {
    kde_peak_week <- fm$kde
    peaks <- fm$x
    x.new <- rnorm(n_sim, 
                   mean=sample(peaks, size = n_sim, replace = TRUE), 
                   sd=kde_peak_week$bw)
    pred_peaks_rounded <- round(x.new)
    
    ## removing early/late peak samples
    idx_out_of_bounds <- which(pred_peaks_rounded<=0 | pred_peaks_rounded>52)
    if(length(idx_out_of_bounds)>0){
      pred_peaks_rounded[-idx_out_of_bounds] 
    }
    
    peaks_binned <- tabulate(pred_peaks_rounded, nbins=52)

    peak_bin_probs <- peaks_binned/(n_sim-length(idx_out_of_bounds))
    names(peak_bin_probs) <- as.character(1:52)
    
    return(peak_bin_probs)    
}


#' Compute predictive distribution for peak week incidence from KDE fit
#'
#' @param fm output from fit_kde_log_peak_week_inc()
#' @param bins cutpoints for incidence bins 
#' @param bin_names vector of bin names, length 1 fewer than length(bins)
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return vector of probabilities representing the predictive distribution of NOT LOG peak week incidence, 
#'         in given bins 
#'
predict_kde_log_peak_week_inc <- function(fm, bins, bin_names, n_sim) {
    kde_peak_week_inc <- fm$kde
    log_peak_inc <- fm$x
    x_new <- rnorm(n_sim, 
                   mean=base::sample(log_peak_inc, size = n_sim, replace = TRUE), 
                   sd=kde_peak_week_inc$bw)
    pred_peaks_binned <- cut(exp(x_new), breaks=bins, right=FALSE) ## CDC incidence targets are [a,b)
    peak_inc_bin_probs <- table(pred_peaks_binned)/n_sim
    
    names(peak_inc_bin_probs) <- bin_names
    return(peak_inc_bin_probs)    
}


#' Compute predictive distribution for weekly incidence from GAM fit
#'
#' @param fm a GAM fit for weekly incidence
#' @param season_weeks season_weeks to predict for
#' @param bins cutpoints for incidence bins 
#' @param bin_names vector of bin names, length 1 fewer than length(bins)
#' @param n_sim number of draws from predictive distribution to return
#'
#' @return matrix of probabilities representing the predictive distribution of NOT LOG scale weekly incidence, 
#'         with row correponding to the bin and column corresponding to each given season_week
#'
predict_kde_log_weekly_inc <- function(fm, season_weeks, bins, bin_names, n_sim) {
    require(mgcv)

    ## code adapted from example at: https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html
    
    ## extract parameter estiamtes and cov matrix...
    beta <- coef(fm)
    Vb <- vcov(fm)
    
    ## simulate replicate beta vectors from posterior...
    Cv <- chol(Vb)
    nb <- length(beta)
    br <- t(Cv) %*% matrix(rnorm(n_sim*nb),nb,n_sim) + beta
    
    ## turn these into replicate linear predictors...
    Xp <- predict(fm,newdata=data.frame(season_week=season_weeks),type="lpmatrix")
    lp <- Xp%*%br
    fv <- lp ## ... finally, replicate expected value vectors
    
    ## now simulate from normal deviates with mean as in fv
    ## and estimated scale...
    tmp <- matrix(rnorm(length(fv), mean=fv, sd=sqrt(fm$sig2)), nrow=nrow(fv), ncol(fv))
    
    #plot(rep(xp,n_sim),exp(tmp),pch=".", ylim=c(-1, 10)) ## plotting replicates
    #points(data$season_week,data$weighted_ili,pch=19,cex=.5) ## and original data
    
    ## compute 95% prediction interval...
    #PI <- apply(exp(tmp),1,quantile,prob=c(.025,0.975))
    
    weekly_inc_bin_probs <- apply(tmp, 
                                  MARGIN=1,
                                  FUN = function(x) {
                                      binned_inc <- cut(exp(x), breaks=bins, right=FALSE)
                                      table(binned_inc)/n_sim
                                  })
    rownames(weekly_inc_bin_probs) <- bin_names
    
    return(weekly_inc_bin_probs)    
}
