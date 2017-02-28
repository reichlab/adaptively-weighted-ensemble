library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(awes)

region_strings <- c("National", paste0("Region", 1:10))

pdf("inst/evaluation/em-stacking/check-em-stacking-test-predictions.pdf", width=10)
for(reg in region_strings) {
  preds <- assemble_predictions(
    preds_path = "inst/evaluation/test-predictions",
    regions = reg,
    models = c("kde", "kcde", "sarima", "em_stacking")
  ) %>%
    mutate(region_season_week =
        paste(region, analysis_time_season, analysis_time_season_week, sep = "_"))
  
  for(prediction_target in c("onset", "peak_week", "peak_inc")) {
    bin_log_prob_cols <- grep(
      paste0(prediction_target, "_", "bin_.*_log_prob"),
      colnames(preds))
    bin_log_prob_cols <- colnames(preds)[bin_log_prob_cols]
    prediction_cols_to_keep <- c("model",
      "region",
      "analysis_time_season",
      "analysis_time_season_week",
      "region_season_week",
      bin_log_prob_cols)
    
    target_preds <- preds %>%
      select_(.dots = prediction_cols_to_keep) %>%
      gather_("bin", "log_prob", bin_log_prob_cols) %>%
      mutate(prob = exp(log_prob),
        bin = as.numeric(substr(bin, 6 + nchar(prediction_target), nchar(bin) - 9)))
    target_preds$bin[is.na(target_preds$bin)] <- 0
    
    temp <- target_preds %>%
      group_by(model, region_season_week) %>%
      summarize(bin_probs_sum = sum(prob))
    if(!isTRUE(all.equal(temp$bin_probs_sum, rep(1, nrow(temp))))) {
      stop("bin probs don't sum to 1")
    }
    
    log_score_col <- paste0(prediction_target, "_log_score")
    prediction_cols_to_keep <- c("model",
      "region",
      "analysis_time_season",
      "analysis_time_season_week",
      "region_season_week",
      log_score_col)
    target_preds_log_score <- preds %>%
      select_(.dots = prediction_cols_to_keep) %>%
      `colnames<-`(c("model",
        "region",
        "analysis_time_season",
        "analysis_time_season_week",
        "region_season_week",
        "log_score")) %>%
      filter(model == "em_stacking") %>%
      mutate(score = exp(log_score))
    
    for(season in unique(target_preds$analysis_time_season)) {
      p <- ggplot() +
        geom_line(aes(x = bin, y = prob, colour = model, linetype = model),
          data = target_preds %>% filter(analysis_time_season == season & analysis_time_season_week %in% 10:40)) +
        geom_hline(aes(yintercept = score, colour = model, linetype = model),
          data = target_preds_log_score %>% filter(analysis_time_season == season & analysis_time_season_week %in% 10:40)) +
        facet_wrap(~ analysis_time_season_week, nrow = 4) +
        ylim(0, 1) +
        ggtitle(paste(reg, season, prediction_target, sep = " - ")) +
        theme_bw()
      print(p)
    }
  }
}
dev.off()


## check summation to 1
for(reg in region_strings) {
  preds <- assemble_predictions(
    preds_path = "inst/evaluation/test-predictions",
    regions = reg,
    models = c("kde", "kcde", "sarima", "em_stacking")
  ) %>%
    mutate(region_season_week =
        paste(region, analysis_time_season, analysis_time_season_week, sep = "_"))
  
  for(prediction_target in c("onset", "peak_week", "peak_inc")) {
    bin_log_prob_cols <- grep(
      paste0(prediction_target, "_", "bin_.*_log_prob"),
      colnames(preds))
    bin_log_prob_cols <- colnames(preds)[bin_log_prob_cols]
    prediction_cols_to_keep <- c("model",
      "region",
      "analysis_time_season",
      "analysis_time_season_week",
      "region_season_week",
      bin_log_prob_cols)
    
    target_preds <- preds %>%
      select_(.dots = prediction_cols_to_keep) %>%
      gather_("bin", "log_prob", bin_log_prob_cols) %>%
      mutate(prob = exp(log_prob),
        bin = as.numeric(substr(bin, 6 + nchar(prediction_target), nchar(bin) - 9)))
    target_preds$bin[is.na(target_preds$bin)] <- 0
    
    temp <- target_preds %>%
      filter(analysis_time_season_week %in% 10:40) %>%
      group_by(model, region_season_week) %>%
      summarize(bin_probs_sum = sum(prob))
    if(!isTRUE(all.equal(temp$bin_probs_sum, rep(1, nrow(temp))))) {
      stop("bin probs don't sum to 1")
    }
  }
}
