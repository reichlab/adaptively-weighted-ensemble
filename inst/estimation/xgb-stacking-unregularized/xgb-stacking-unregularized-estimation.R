### Estimate function that gives model weights based on observed inputs

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xgboost)
library(xgbstack)
library(awes)

loso_preds_path <- "inst/estimation/loso-predictions"
stacking_model_fits_path <- "inst/estimation/xgb-stacking-unregularized/fits"
component_model_names <- c("kde", "kcde", "sarima")

### Command line arguments
args <- commandArgs(trailingOnly = TRUE)
## "National" or "Regionk" for k in {1, ..., 10}
region <- args[1]

## "onset", "peak_week", "peak_inc"
prediction_target <- args[2]

## Subset of the following, separated by "-"
## "analysis_time_season_week" "kcde_model_confidence" "sarima_model_confidence" "weighted_ili"
## explanatory_variables <- c("analysis_time_season_week", "kcde_model_confidence", "sarima_model_confidence")
explanatory_variables <- strsplit(args[3], "-")[[1]]

## integer number of threads to use
nthread <- as.integer(args[4])

### get xgbstack fit
xgbstack_fit <-
  fit_stacked_model(
    regions = region,
    prediction_target = prediction_target,
    component_model_names = component_model_names,
    explanatory_variables = explanatory_variables,
    loso_preds_path = loso_preds_path,
    min_child_weight = - 10^10,
    eta = 0.001,
    max_depth = 10,
    lambda = 0,
    alpha = 0,
    gamma = 0,
    nrounds = 100000,
    verbose = 1,
    nthread = nthread
  )

## get resulting model weights at a grid of values for chosen explanatory variables
if("kcde_model_confidence" %in% explanatory_variables ||
    "sarima_model_confidence" %in% explanatory_variables) {
  if(identical(prediction_target, "peak_inc")) {
    model_confidence_grid <- seq(from = 1, to = 99, by = 2)
  } else if(identical(prediction_target, "onset")) {
    model_confidence_grid <- seq_len(34)
  } else {
    ## peak_week
    model_confidence_grid <- seq_len(33)
  }
} else {
  model_confidence_grid <- 1
}

if("weighted_ili" %in% explanatory_variables) {
  weighted_ili_grid <- seq(from = 0, to = 13, by = 0.5)
} else {
  weighted_ili_grid <- 1
}

typical_analysis_time_season_week <- 17L
if(prediction_target %in% c("onset", "peak_week")) {
  typical_model_confidence <- 5L
} else {
  typical_model_confidence <- 20L
}
typical_ili <- 2.5

newdata <- bind_rows(
  expand.grid(
    analysis_time_season_week = 10:40,
    kcde_model_confidence = model_confidence_grid,
    sarima_model_confidence = typical_model_confidence,
    weighted_ili = typical_ili,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  expand.grid(
    analysis_time_season_week = 10:40,
    kcde_model_confidence = typical_model_confidence,
    sarima_model_confidence = model_confidence_grid,
    weighted_ili = typical_ili,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  expand.grid(
    analysis_time_season_week = 10:40,
    kcde_model_confidence = typical_model_confidence,
    sarima_model_confidence = typical_model_confidence,
    weighted_ili = weighted_ili_grid,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  expand.grid(
    analysis_time_season_week = typical_analysis_time_season_week,
    kcde_model_confidence = model_confidence_grid,
    sarima_model_confidence = model_confidence_grid,
    weighted_ili = typical_ili,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  expand.grid(
    analysis_time_season_week = typical_analysis_time_season_week,
    kcde_model_confidence = model_confidence_grid,
    sarima_model_confidence = typical_model_confidence,
    weighted_ili = weighted_ili_grid,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  expand.grid(
    analysis_time_season_week = typical_analysis_time_season_week,
    kcde_model_confidence = typical_model_confidence,
    sarima_model_confidence = model_confidence_grid,
    weighted_ili = weighted_ili_grid,
    stringsAsFactors = FALSE)[, explanatory_variables, drop = FALSE],
  assemble_stacking_inputs(
    regions = c("National", paste0("Region", 1:10)),
    prediction_target = prediction_target,
    component_model_names = component_model_names,
    explanatory_variables = explanatory_variables,
    include_model_performance = FALSE,
    preds_path = "inst/estimation/loso-predictions"
  )[, explanatory_variables, drop = FALSE],
  assemble_stacking_inputs(
    regions = c("National", paste0("Region", 1:10)),
    prediction_target = prediction_target,
    component_model_names = component_model_names,
    explanatory_variables = explanatory_variables,
    include_model_performance = FALSE,
    preds_path = "inst/evaluation/test-predictions"
  )[, explanatory_variables, drop = FALSE]
) %>%
  distinct()


## save results over a range of values of nrounds; used to ensure convergence
ntreelimit_vals = c(50000, 60000, 70000, 80000, 90000, 95000, 97000, 98000, 99000, 100000)

for(ntreelimit in ntreelimit_vals) {
  model_weights <- compute_model_weights(
    xgbstack_fit,
    newdata = newdata,
    ntreelimit = ntreelimit,
    log = FALSE,
    format = "bare")
  model_weights <- cbind(newdata, model_weights)
  
  saveRDS(model_weights,
    file = file.path(stacking_model_fits_path,
      paste0("model_weights_",
        region, "_",
        prediction_target, "_",
        paste0(explanatory_variables, collapse = "-"), "_",
        ntreelimit,
        ".rds")))
}
