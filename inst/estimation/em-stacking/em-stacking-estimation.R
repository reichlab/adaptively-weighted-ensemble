### Estimate function that gives model weights based on observed inputs

library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
library(awes)

registerDoMC(4)

models <- c("kde", "kcde", "sarima")

## estimate model weights for each prediction target, region, and analysis_time_season_week
weights <- foreach(region = c("National", paste0("Region", 1:10))) %dopar% {
  rbind.fill(lapply(
    c("onset", "peak_week", "peak_inc"),
    function(prediction_target) {
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
      
      ## get weights combining across all analysis_time_season_weeks
      const_model_weights <- fit_unregularized_stacked_model(
        component_model_log_scores = loso_pred_res[, models],
        method = "em",
        tol = .Machine$double.eps) %>%
        matrix(nrow = 1) %>%
        `colnames<-`(models) %>%
        as.data.frame() %>%
        mutate(region = region,
          prediction_target = prediction_target,
          analysis_time_season_week = "all-combined")
      
      ### separate model weights by season week, combining across all seasons
      weekly_model_weights <- sapply(
        unique(loso_pred_res$analysis_time_season_week),
        function(season_week) {
          fit_unregularized_stacked_model(
            component_model_log_scores = loso_pred_res[
              loso_pred_res$analysis_time_season_week == season_week, models],
            method = "em",
            tol = .Machine$double.eps)
        }) %>%
        t() %>%
        `colnames<-`(models) %>%
        as.data.frame() %>%
        mutate(
          region = region,
          prediction_target = prediction_target,
          analysis_time_season_week = as.character(unique(loso_pred_res$analysis_time_season_week))
        )
      
      return(bind_rows(const_model_weights, weekly_model_weights))
    }))
}
weights <- rbind.fill(weights)
weights <- weights[, c("region", "prediction_target", "analysis_time_season_week", "kde", "kcde", "sarima")]

saveRDS(weights, "inst/estimation/em-stacking/fits/em_weights.rds")
