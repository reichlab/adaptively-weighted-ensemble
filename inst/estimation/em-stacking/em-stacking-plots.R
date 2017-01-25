### Estimate function that gives model weights based on observed inputs

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(xgboost)
library(xgbstack)
library(doMC)
library(awes)

pdf("inst/estimation/em-stacking/plots/em-stacking-plots-all-regions-all-targets.pdf")

weights <- readRDS("inst/estimation/em-stacking/fits/em_weights.rds") %>%
  gather_("model", "weight", c("kde", "kcde", "sarima"))

models <- c("kde", "kcde", "sarima")
for(region in c("National", paste0("Region", 1:10))) {
  for(prediction_target in c("onset", "peak_week", "peak_inc")) {
    loso_pred_res <- assemble_loso_predictions(
      loso_preds_path = "inst/estimation/loso-predictions",
      regions = region,
      models = models,
      prediction_targets = prediction_target,
      prediction_types = "log_score"
    ) %>%
      spread_("model", paste0(prediction_target, "_log_score"))
    
    ## drop rows where any of the methods have NA values.  this is agressive
    loso_pred_res <- loso_pred_res[
      apply(loso_pred_res[, models],
        1,
        function(x) {!any(is.na(x))}), # row has na?  drop if so
      , # all columns
      drop = FALSE]
    
    mean_loso_pred_res <- loso_pred_res %>%
      mutate(
        kcde = sapply(kcde, function(x) max(x, -10)),
        kde = sapply(kde, function(x) max(x, -10)),
        sarima = sapply(sarima, function(x) max(x, -10))) %>%
      gather("model", "log_score", kcde, kde, sarima) %>%
      group_by(model, analysis_time_season_week) %>%
      summarize(mean_log_score = mean(log_score)) %>%
      as.data.frame()
    
    mean_log_scores_plot <- ggplot(mean_loso_pred_res) +
      geom_line(aes(x = analysis_time_season_week, y = mean_log_score, colour = model)) +
      theme_bw()
    
    model_weights_plot <- ggplot() +
      geom_line(aes(x = as.integer(analysis_time_season_week), y = weight, colour = model),
        data = weights[weights$region == region &
            weights$prediction_target == prediction_target &
            weights$analysis_time_season_week != "all-combined", ]) +
      geom_hline(aes(yintercept = weight, colour = model),
        linetype = 2,
        weights[weights$region == region &
            weights$prediction_target == prediction_target &
            weights$analysis_time_season_week == "all-combined", ]) +
      theme_bw()
    
    grid.newpage()
    pushViewport(viewport(layout =
        grid.layout(nrow = 3, ncol = 1,
          heights = unit(rep(1, 3), c("lines", "null", "null")))))
    grid.text(paste0(region, " - ", prediction_target), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(mean_log_scores_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(model_weights_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  }
}

dev.off()
