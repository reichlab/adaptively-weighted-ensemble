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

source("R/stacking-utils.R")

loso_preds_path <- "inst/estimation/loso-predictions"
stacking_model_fits_path <- "inst/estimation/stacking/stacking-fits-dev"

region <- "National"
# prediction_targets <- c("onset",
#   "peak_week",
#   "peak_inc",
#   "ph_1_inc",
#   "ph_2_inc",
#   "ph_3_inc",
#   "ph_4_inc")
prediction_target <- "peak_week"

component_model_names <- c("kde", "kcde", "sarima")

## Assemble training data
regions <- region
loso_pred_res <- rbind.fill(lapply(
  Sys.glob(paste0(loso_preds_path, "/", component_model_names, "-", regions, "*")),
  readRDS
))

target_loso_pred_res <- loso_pred_res %>%
  select_("model",
    "analysis_time_season",
    "analysis_time_season_week",
    paste0(prediction_target, "_log_score")) %>%
  spread_("model", paste0(prediction_target, "_log_score"))

## drop rows where any of the methods have NA values.  this is agressive
target_loso_pred_res <- target_loso_pred_res[
  apply(target_loso_pred_res[, component_model_names],
    1,
    function(x) {!any(is.na(x))}), # row has na?  drop if so
  , # all columns
  drop = FALSE]

mean_target_loso_pred_res <- target_loso_pred_res %>%
  mutate(
    kcde = sapply(kcde, function(x) max(x, -10)),
    kde = sapply(kde, function(x) max(x, -10)),
    sarima = sapply(sarima, function(x) max(x, -10))) %>%
  gather("model", "log_score", kcde, kde, sarima) %>%
  group_by(model, analysis_time_season_week) %>%
  summarize(mean_log_score = mean(log_score)) %>%
  as.data.frame()


### single model weights across all analysis_time_season_week values...
## ...and all seasons
single_model_weights <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "em",
  tol = .Machine$double.eps) %>%
  `names<-`(c("kcde", "kde", "sarima"))

single_model_weights_safe <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "em-safe",
  tol = .Machine$double.eps) %>%
  `names<-`(c("kcde", "kde", "sarima"))

single_model_weights_grid <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "grid-search") %>%
  `names<-`(c("kcde", "kde", "sarima"))

## ...separately for each training season
single_model_weights_by_season <- sapply(
  unique(target_loso_pred_res$analysis_time_season),
  function(season) {
    fit_unregularized_stacked_model(
    component_model_log_scores = target_loso_pred_res[
      target_loso_pred_res$analysis_time_season == season, c("kcde", "kde", "sarima")],
    method = "em",
    tol = .Machine$double.eps)
  }) %>%
  t() %>%
  `colnames<-`(c("kcde", "kde", "sarima"))

### separate model weights by season week, combining across all seasons
week_model_weights <- sapply(
  unique(target_loso_pred_res$analysis_time_season_week),
  function(season_week) {
    fit_unregularized_stacked_model(
      component_model_log_scores = target_loso_pred_res[
        target_loso_pred_res$analysis_time_season_week == season_week, c("kcde", "kde", "sarima")],
      method = "em",
      tol = .Machine$double.eps)
  }) %>%
  t() %>%
  `colnames<-`(c("kcde", "kde", "sarima")) %>%
  as.data.frame() %>%
  mutate(
    analysis_time_season_week = 10:40
  ) %>%
  gather("model", "weight", kcde, kde, sarima)

### separate model weights by season week, separately for each season
week_model_weights_by_season <- rbind.fill(lapply(unique(target_loso_pred_res$analysis_time_season),
  function(season) {
    sapply(
      unique(target_loso_pred_res$analysis_time_season_week[target_loso_pred_res$analysis_time_season == season]),
      function(season_week) {
        fit_unregularized_stacked_model(
          component_model_log_scores = target_loso_pred_res[
            target_loso_pred_res$analysis_time_season == season &
              target_loso_pred_res$analysis_time_season_week == season_week,
            c("kcde", "kde", "sarima")],
          method = "em",
          tol = .Machine$double.eps)
      }) %>%
      t() %>%
      `colnames<-`(c("kcde", "kde", "sarima")) %>%
      as.data.frame() %>%
      mutate(
        analysis_time_season_week = 
          unique(target_loso_pred_res$analysis_time_season_week[target_loso_pred_res$analysis_time_season == season])
      ) %>%
      gather("model", "weight", kcde, kde, sarima) %>%
      mutate(analysis_time_season = season)
  }))



pdf("inst/estimation/stacking/regularization-exploration/unregularized_stacking_plots_National_peak_week.pdf")

ggplot(mean_target_loso_pred_res) +
  geom_line(aes(x = analysis_time_season_week, y = mean_log_score, colour = model)) +
  xlab("Week of Season at Analysis Time") +
  ylab("Mean Cross-Validated Log Score") +
  ggtitle("Mean Log Scores by Model and Season Week") +
  theme_bw()

combined_single_model_weights <- bind_rows(
  single_model_weights_by_season %>%
    as.data.frame() %>%
    mutate(season = rownames(single_model_weights_by_season)) %>%
    gather("model", "weight", seq_len(ncol(single_model_weights_by_season))),
  single_model_weights %>%
    as.matrix() %>%
    as.data.frame() %>%
    `colnames<-`("weight") %>%
    mutate(model = names(single_model_weights),
      season = "all seasons")
)
  
ggplot() +
  geom_point(aes(x = season, y = weight, colour = model),
    data = combined_single_model_weights) +
  geom_hline(aes(yintercept = weight, colour = model),
    data = combined_single_model_weights[combined_single_model_weights$season == "all seasons", ]) +
  theme_bw()

ggplot() +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    data = week_model_weights) +
  geom_hline(aes(yintercept = weight, colour = model),
    linetype = 2,
    data = single_model_weights %>%
      as.matrix() %>%
      as.data.frame() %>%
      `colnames<-`("weight") %>%
      mutate(model = names(single_model_weights))) +
  xlab("Week of Season at Analysis Time") +
  ylab("Component Model Weights") +
  ggtitle("Unregularized Model Weights by Season Week") +
  theme_bw()

ggplot() +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    alpha = 0.2,
    data = week_model_weights_by_season) +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    data = week_model_weights) +
  facet_wrap( ~ analysis_time_season) +
  xlab("Week of Season at Analysis Time") +
  ylab("Component Model Weights") +
  ggtitle("Unregularized Model Weights by Season Week") +
  theme_bw()

dev.off()






region <- c("National", paste0("Region", 1:10))
# prediction_targets <- c("onset",
#   "peak_week",
#   "peak_inc",
#   "ph_1_inc",
#   "ph_2_inc",
#   "ph_3_inc",
#   "ph_4_inc")
prediction_target <- "peak_week"

component_model_names <- c("kde", "kcde", "sarima")

## Assemble training data
regions <- region
loso_pred_res <- rbind.fill(lapply(
  Sys.glob(outer(
    component_model_names,
    regions,
    function(model, region) {
      paste0(loso_preds_path, "/", model, "-", region, "-*")
    }
  )),
  function(file_path) {
    region_val <- names(unlist(sapply(regions, function(region_val) grep(paste0(region_val, "-"), file_path))))
    readRDS(file_path) %>%
      mutate(region = region_val)
  }
))

target_loso_pred_res <- loso_pred_res %>%
  select_("model",
    "region",
    "analysis_time_season",
    "analysis_time_season_week",
    paste0(prediction_target, "_log_score")) %>%
  spread_("model", paste0(prediction_target, "_log_score"))

## drop rows where any of the methods have NA values.  this is agressive
target_loso_pred_res <- target_loso_pred_res[
  apply(target_loso_pred_res[, component_model_names],
    1,
    function(x) {!any(is.na(x))}), # row has na?  drop if so
  , # all columns
  drop = FALSE]

mean_target_loso_pred_res <- target_loso_pred_res %>%
  mutate(
    kcde = sapply(kcde, function(x) max(x, -10)),
    kde = sapply(kde, function(x) max(x, -10)),
    sarima = sapply(sarima, function(x) max(x, -10))) %>%
  gather("model", "log_score", kcde, kde, sarima) %>%
  group_by(model, analysis_time_season_week) %>%
  summarize(mean_log_score = mean(log_score)) %>%
  as.data.frame()


### single model weights across all analysis_time_season_week values...
## ...and all seasons
single_model_weights <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "em",
  tol = .Machine$double.eps) %>%
  `names<-`(c("kcde", "kde", "sarima"))

single_model_weights_safe <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "em-safe",
  tol = .Machine$double.eps) %>%
  `names<-`(c("kcde", "kde", "sarima"))

single_model_weights_grid <- fit_unregularized_stacked_model(
  component_model_log_scores = target_loso_pred_res[, c("kcde", "kde", "sarima")],
  method = "grid-search") %>%
  `names<-`(c("kcde", "kde", "sarima"))

## ...separately for each training season
single_model_weights_by_season <- sapply(
  unique(target_loso_pred_res$analysis_time_season),
  function(season) {
    fit_unregularized_stacked_model(
      component_model_log_scores = target_loso_pred_res[
        target_loso_pred_res$analysis_time_season == season, c("kcde", "kde", "sarima")],
      method = "em",
      tol = .Machine$double.eps)
  }) %>%
  t() %>%
  `colnames<-`(c("kcde", "kde", "sarima"))

### separate model weights by season week, combining across all seasons
week_model_weights <- sapply(
  unique(target_loso_pred_res$analysis_time_season_week),
  function(season_week) {
    fit_unregularized_stacked_model(
      component_model_log_scores = target_loso_pred_res[
        target_loso_pred_res$analysis_time_season_week == season_week, c("kcde", "kde", "sarima")],
      method = "em",
      tol = .Machine$double.eps)
  }) %>%
  t() %>%
  `colnames<-`(c("kcde", "kde", "sarima")) %>%
  as.data.frame() %>%
  mutate(
    analysis_time_season_week = 10:40
  ) %>%
  gather("model", "weight", kcde, kde, sarima)

### separate model weights by season week, separately for each season
week_model_weights_by_season <- rbind.fill(lapply(unique(target_loso_pred_res$analysis_time_season),
  function(season) {
    sapply(
      unique(target_loso_pred_res$analysis_time_season_week[target_loso_pred_res$analysis_time_season == season]),
      function(season_week) {
        fit_unregularized_stacked_model(
          component_model_log_scores = target_loso_pred_res[
            target_loso_pred_res$analysis_time_season == season &
              target_loso_pred_res$analysis_time_season_week == season_week,
            c("kcde", "kde", "sarima")],
          method = "em",
          tol = .Machine$double.eps)
      }) %>%
      t() %>%
      `colnames<-`(c("kcde", "kde", "sarima")) %>%
      as.data.frame() %>%
      mutate(
        analysis_time_season_week = 
          unique(target_loso_pred_res$analysis_time_season_week[target_loso_pred_res$analysis_time_season == season])
      ) %>%
      gather("model", "weight", kcde, kde, sarima) %>%
      mutate(analysis_time_season = season)
  }))



pdf("inst/estimation/stacking/regularization-exploration/unregularized_stacking_plots_all_regions_peak_week.pdf")

ggplot(mean_target_loso_pred_res) +
  geom_line(aes(x = analysis_time_season_week, y = mean_log_score, colour = model)) +
  xlab("Week of Season at Analysis Time") +
  ylab("Mean Cross-Validated Log Score") +
  ggtitle("Mean Log Scores by Model and Season Week") +
  theme_bw()

combined_single_model_weights <- bind_rows(
  single_model_weights_by_season %>%
    as.data.frame() %>%
    mutate(season = rownames(single_model_weights_by_season)) %>%
    gather("model", "weight", seq_len(ncol(single_model_weights_by_season))),
  single_model_weights %>%
    as.matrix() %>%
    as.data.frame() %>%
    `colnames<-`("weight") %>%
    mutate(model = names(single_model_weights),
      season = "all seasons")
)

ggplot() +
  geom_point(aes(x = season, y = weight, colour = model),
    data = combined_single_model_weights) +
  geom_hline(aes(yintercept = weight, colour = model),
    data = combined_single_model_weights[combined_single_model_weights$season == "all seasons", ]) +
  theme_bw()

ggplot() +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    data = week_model_weights) +
  geom_hline(aes(yintercept = weight, colour = model),
    linetype = 2,
    data = single_model_weights %>%
      as.matrix() %>%
      as.data.frame() %>%
      `colnames<-`("weight") %>%
      mutate(model = names(single_model_weights))) +
  xlab("Week of Season at Analysis Time") +
  ylab("Component Model Weights") +
  ggtitle("Unregularized Model Weights by Season Week") +
  theme_bw()

ggplot() +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    data = week_model_weights_by_season) +
  geom_line(aes(x = analysis_time_season_week, y = weight, colour = model),
    alpha = 0.4,
    data = week_model_weights) +
  facet_wrap( ~ analysis_time_season) +
  xlab("Week of Season at Analysis Time") +
  ylab("Component Model Weights") +
  ggtitle("Unregularized Model Weights by Season Week") +
  theme_bw()

dev.off()







pdf("inst/estimation/stacking/regularization-exploration/unregularized_stacking_plots_all_regions_all_targets.pdf")

weights <- readRDS("inst/estimation/stacking/fits-unregularized/unregularized_weights.rds") %>%
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
