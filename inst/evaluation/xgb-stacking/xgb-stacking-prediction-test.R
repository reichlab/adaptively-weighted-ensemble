### Estimate function that gives model weights based on observed inputs

library(plyr)
library(dplyr)
library(tidyr)
library(xgboost)
library(xgbstack)
library(doMC)
library(awes)

loso_preds_path <- "inst/estimation/loso-predictions"
stacking_model_fits_path <- "inst/estimation/xgb-stacking/fits"

component_model_names <- c("kde", "kcde", "sarima")

all_explanatory_variables_combos <- c(
  "analysis_time_season_week",
  "analysis_time_season_week-kcde_model_confidence-sarima_model_confidence",
  "analysis_time_season_week-kcde_model_confidence-sarima_model_confidence-weighted_ili"
)

registerDoMC(4)

foreach(region = c("National", paste0("Region", 1:10))) %dopar% {
  # region <- "National"
  # region <- "Region6"
  for(explanatory_variables in all_explanatory_variables_combos) {
    explanatory_variables_split <- strsplit(explanatory_variables, "-")[[1]]
    if(length(explanatory_variables_split) == 1) {
      ensemble_model_name <- "xgb_stacking_reg_w"
    } else if(length(explanatory_variables_split) == 3) {
      ensemble_model_name <- "xgb_stacking_reg_wu"
    } else {
      ensemble_model_name <- "xgb_stacking_reg_wui"
    }
    
    weights <- NULL
    for(target in c("onset", "peak_week", "peak_inc")) {
      target_weights <- readRDS(
        file = file.path(stacking_model_fits_path,
          paste0("model_weights_", region, "_", target, "_", explanatory_variables, ".rds"))) %>%
          select_(.dots = c(explanatory_variables_split, paste0(component_model_names, "_log_score_params_combined"))) %>%
        `colnames<-`(c(explanatory_variables_split, component_model_names)) %>%
        mutate(region = region,
          prediction_target = target)
      
      weights <- bind_rows(weights, target_weights)
    }
    
    weights <- weights %>%
      gather_("model", "weight", c("kde", "kcde", "sarima")) %>%
      spread("prediction_target", "weight")
    
    for(season in paste0(2011:2015, "/", 2012:2016)) {
      # season <- "2011/2012"
      ensemble_predictions <- NULL
      for(prediction_target in c("onset", "peak_week", "peak_inc")) {
        # prediction_target <- "peak_week"
        stacking_features <-
          assemble_stacking_inputs(
            regions = region,
            seasons = season,
            prediction_target = prediction_target,
            component_model_names = component_model_names,
            explanatory_variables = explanatory_variables_split,
            include_model_performance = FALSE,
            preds_path = "inst/evaluation/test-predictions"
          ) %>%
          select_(
            .dots = c("analysis_time_season", "analysis_time_season_week", explanatory_variables_split)
          )
        component_preds <- assemble_predictions(
          preds_path = "inst/evaluation/test-predictions",
          regions = region,
          seasons = season,
          prediction_targets = prediction_target,
          models = component_model_names
        )
        component_preds <- left_join(
          component_preds,
          stacking_features,
          by = c("analysis_time_season", "analysis_time_season_week")
        )
        
        ensemble_predictions_target <- weighted_combine_predictions(
          component_preds = component_preds,
          weights = weights)
        if(is.null(ensemble_predictions)) {
          ensemble_predictions <- ensemble_predictions_target
        } else {
          ensemble_predictions <- left_join(
            ensemble_predictions,
            ensemble_predictions_target,
            by = c("region", "analysis_time_season", "analysis_time_season_week")
          )
        }
      }
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
}
