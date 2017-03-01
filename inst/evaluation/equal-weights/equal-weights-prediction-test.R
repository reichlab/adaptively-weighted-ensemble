library(plyr)
library(dplyr)
library(tidyr)
library(doMC)
library(awes)

registerDoMC(3)

models <- c("kde", "kcde", "sarima")
ensemble_model_name <- "equal_weights"

## data frame of weights with columns named region, model, onset, peak_inc, peak_week
weights <- data.frame(
  model = models,
  onset = 1 / length(models),
  peak_inc = 1 / length(models),
  peak_week = 1 / length(models)
)

## estimate model weights for each prediction target, region, and analysis_time_season_week
foreach(region = c("National", paste0("Region", 1:10))) %dopar% {
  # region <- "National"
  # region <- "Region1"
  for(season in paste0(2011:2015, "/", 2012:2016)) {
    # season <- "2011/2012"
    component_preds <- assemble_predictions(
      preds_path = "inst/evaluation/test-predictions",
      regions = region,
      seasons = season,
      models = models
    )
    
    ensemble_predictions <- weighted_combine_predictions(
      component_preds = component_preds,
      weights = weights)
    ensemble_predictions <- cbind(
      data.frame(model = ensemble_model_name),
      ensemble_predictions
    )
    
    saveRDS(ensemble_predictions,
      paste0("inst/evaluation/test-predictions/",
        ensemble_model_name, "-",
        region, "-",
        gsub("/", "-", season),
        "-test-predictions.rds")
    )
  }
}
