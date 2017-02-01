library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(awes)

region_strings <- c("National", paste0("Region", 1:10))

pdf("inst/estimation/sarima/check-sarima-predictions.pdf", width=10)
for(reg in region_strings) {
  loso_preds <- assemble_loso_predictions(
    regions = reg,
    models = "sarima")
  
  temp <- loso_preds %>%
    select(model,
      region,
      analysis_time_season,
      analysis_time_season_week,
      onset_log_score,
      peak_week_log_score,
      peak_inc_log_score,
      ph_1_inc_log_score,
      ph_2_inc_log_score,
      ph_3_inc_log_score,
      ph_4_inc_log_score
    ) %>%
    gather("prediction_target", "log_score", 5:11) %>%
    mutate(
      prediction_target = substr(prediction_target, 1, nchar(prediction_target) - 10))
  
  p <- ggplot(temp,
    aes(x = analysis_time_season_week, y = log_score)) +
    geom_line(aes(color = analysis_time_season)) +
    facet_grid(. ~ prediction_target) +
    geom_smooth(se = FALSE, color = "black") +
    ylim(-10, 0) +
    ggtitle(paste0(reg, " - log score"))
  print(p)
  
  temp <- loso_preds %>%
    select(model,
      region,
      analysis_time_season,
      analysis_time_season_week,
      onset_competition_log_score,
      peak_week_competition_log_score,
      peak_inc_competition_log_score,
      ph_1_inc_competition_log_score,
      ph_2_inc_competition_log_score,
      ph_3_inc_competition_log_score,
      ph_4_inc_competition_log_score
    ) %>%
    gather("prediction_target", "log_score", 5:11) %>%
    mutate(
      prediction_target = substr(prediction_target, 1, nchar(prediction_target) - 10))
  
  p <- ggplot(temp,
    aes(x = analysis_time_season_week, y = log_score)) +
    geom_line(aes(color = analysis_time_season)) +
    facet_grid(. ~ prediction_target) +
    geom_smooth(se = FALSE, color = "black") +
    ylim(-10, 0) +
    ggtitle(paste0(reg, " - competition log score"))
  print(p)
}
dev.off()


for(prediction_target in c("onset", "peak_week", "peak_inc")) {
  pdf(paste0("inst/estimation/sarima/check-sarima-predictions-", prediction_target, ".pdf"), width=10)
  for(reg in region_strings) {
    loso_preds <- assemble_loso_predictions(
      regions = reg,
      prediction_target = prediction_target,
      models = "sarima")
    
    cols_to_examine <- grep(paste0(prediction_target, ".*_log_prob"), colnames(loso_preds), value = TRUE)
    
    temp <- loso_preds %>%
      select_(.dots = c("model",
        "region",
        "analysis_time_season",
        "analysis_time_season_week",
        cols_to_examine)) %>%
      gather_("bin", "log_prob", cols_to_examine) %>%
      mutate(bin = as.numeric(substr(bin, nchar(prediction_target) + 6, nchar(bin) - 9)),
        prob = exp(log_prob))
    if(identical(prediction_target, "onset")) {
      temp$bin[is.na(temp$bin)] <- 0
    }
    
    ## Load and clean up data set
    data <- read.csv("data-raw/allflu-cleaned.csv", stringsAsFactors = FALSE)
    data$time <- as.POSIXct(data$time)
    
    ## subset data to be only the region of interest
    reg_data <- if(reg == "National") {
        "X"
      } else {
        paste0("Region ", substr(reg, 7, nchar(reg)))
      }
    data <- data[data$region == reg_data,]
    
    ## get seasonal quantities
    observed_seasonal_quantities_by_season <- t(sapply(
      paste0(1997:2010, "/", 1998:2011),
      function(season) {
        get_observed_seasonal_quantities(
          data = data,
          season = season,
          first_CDC_season_week = 10,
          last_CDC_season_week = 41,
          onset_baseline = 
            get_onset_baseline(region = reg, season = season),
          incidence_var = "weighted_ili",
          incidence_bins = data.frame(
            lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
            upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
          incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
        )
      }
    ))
    
    if(identical(prediction_target, "onset")) {
      obs_target_by_season <- data.frame(
        analysis_time_season = rep(rownames(observed_seasonal_quantities_by_season),
          times = sapply(observed_seasonal_quantities_by_season[, which(colnames(observed_seasonal_quantities_by_season) == "observed_onset_week")], length)),
        bin = as.numeric(unlist(observed_seasonal_quantities_by_season[, which(colnames(observed_seasonal_quantities_by_season) == "observed_onset_week")]))
      )
      obs_target_by_season$bin[is.na(obs_target_by_season$bin)] <- 0
    } else if(identical(prediction_target, "peak_week")) {
      obs_target_by_season <- data.frame(
        analysis_time_season = rep(rownames(observed_seasonal_quantities_by_season),
          times = sapply(observed_seasonal_quantities_by_season[, which(colnames(observed_seasonal_quantities_by_season) == "observed_peak_week")], length)),
        bin = as.numeric(unlist(observed_seasonal_quantities_by_season[, which(colnames(observed_seasonal_quantities_by_season) == "observed_peak_week")]))
      )
    } else if(identical(prediction_target, "peak_inc")) {
      obs_target_by_season <- data.frame(
        analysis_time_season = rownames(observed_seasonal_quantities_by_season),
        bin = as.numeric(observed_seasonal_quantities_by_season[, which(colnames(observed_seasonal_quantities_by_season) == "observed_peak_inc_bin")])
      )
    }
    
    # loso_preds <- right_join(loso_preds,
    #   obs_target_by_season,
    #   by = c("region",
    #     "analysis_time_season",
    #     "analysis_time_season_week"))
    
    p <- ggplot(temp) +
      geom_line(aes(x = bin, y = prob, alpha = analysis_time_season_week, group = analysis_time_season_week)) +
      geom_vline(aes(xintercept = bin), colour = "red", linetype = 2, data = obs_target_by_season) +
      facet_wrap( ~ analysis_time_season) +
#      geom_smooth(se = FALSE, color = "black") +
#      ylim(-10, 0) +
      ggtitle(paste0(reg, " - ", prediction_target, ifelse(identical(prediction_target, "onset"), " - bin 0 represents 'no onset'", ""))) +
      theme_bw()
    print(p)
  }
  dev.off()
}

for(prediction_target in c("ph_1_inc", "ph_2_inc", "ph_3_inc", "ph_4_inc")) {
  ph <- as.numeric(substr(prediction_target, 4, 4))
  pdf(paste0("inst/estimation/sarima/check-sarima-predictions-", prediction_target, ".pdf"), width=10)
  for(reg in region_strings) {
    loso_preds <- assemble_loso_predictions(
      regions = reg,
      prediction_target = prediction_target,
      models = "sarima")
    
    for(analysis_time_season_to_plot in unique(loso_preds$analysis_time_season)) {
      cols_to_examine <- grep(paste0(prediction_target, ".*_log_prob"), colnames(loso_preds), value = TRUE)
      
      temp <- loso_preds %>%
        select_(.dots = c("model",
          "region",
          "analysis_time_season",
          "analysis_time_season_week",
          cols_to_examine)) %>%
        filter(analysis_time_season == analysis_time_season_to_plot) %>%
        gather_("bin", "log_prob", cols_to_examine) %>%
        mutate(bin = as.numeric(substr(bin, nchar(prediction_target) + 6, nchar(bin) - 9)),
          prob = exp(log_prob))
      
      ## Load and clean up data set
      data <- read.csv("data-raw/allflu-cleaned.csv", stringsAsFactors = FALSE)
      data$time <- as.POSIXct(data$time)
      
      ## subset data to be only the region of interest
      reg_data <- if(reg == "National") {
        "X"
      } else {
        paste0("Region ", substr(reg, 7, nchar(reg)))
      }
      data <- data[data$region == reg_data,]
      data$analysis_time_season <- data$season
      data$analysis_time_season_week <- data$season_week
      data$lead_weighted_ili <- lead(data$weighted_ili, ph)
      
      p <- ggplot(temp) +
        geom_line(aes(x = bin, y = prob)) +
        geom_vline(aes(xintercept = lead_weighted_ili), colour = "red", linetype = 2,
          data = data[data$season == analysis_time_season_to_plot & data$analysis_time_season_week %in% 10:41, ]) +
        facet_wrap( ~ analysis_time_season_week) +
        #      geom_smooth(se = FALSE, color = "black") +
        #      ylim(-10, 0) +
        ggtitle(paste0(reg, " - ", prediction_target, "season = ", analysis_time_season_to_plot, "; faceted by analysis time season week")) +
        theme_bw()
      print(p)
    }
  }
  dev.off()
}
