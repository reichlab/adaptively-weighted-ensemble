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
