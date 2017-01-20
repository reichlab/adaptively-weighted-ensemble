#' Assemble leave-one-season-out predictions made by kde, kcde, and sarima
#' models on training data.
#' 
#' @param loso_preds_path path to leave-one-season-out predictions.
#' @param regions character vector specifying regions for which to get predictions,
#'   "National" or "Regionk" for k in 1, ..., 10.
#' @param models character vector specifying models for which to get predictions,
#'   "kde", "kcde", or "sarima"
#' @param prediction_targets character vector specifying prediction targets,
#'   "onset", "peak_week", "peak_inc", "ph_1_inc", ..., "ph_4_inc"
#' @param prediction_types character vector specifying prediction types,
#'   "log_score", "competition_log_score", or "bin_log_probs"
#' 
#' @return a data frame with predictions.
#' 
#' @export
assemble_loso_predictions <- function(
  loso_preds_path = "inst/estimation/loso-predictions",
  regions = c("National", paste0("Region", 1:10)),
  models = c("kde", "kcde", "sarima"),
  prediction_targets = c("onset", "peak_week", "peak_inc", "ph_1_inc", "ph_2_inc", "ph_3_inc", "ph_4_inc"),
  prediction_types = c("log_score", "competition_log_score", "bin_log_probs")
  ) {
  
  loso_pred_res <- rbind.fill(lapply(
    Sys.glob(outer(
      models,
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
  
  prediction_cols_to_keep_templates <- outer(prediction_targets, prediction_types,
    function(target, type) {
      paste0(target, "_", ifelse(type == "bin_log_probs", "bin_.*_log_prob", type))
    }) %>%
    as.vector()
  prediction_cols_to_keep <- lapply(
    prediction_cols_to_keep_templates,
    function(pattern) grep(pattern, names(loso_pred_res))) %>%
    unlist()
  prediction_cols_to_keep <- names(loso_pred_res)[prediction_cols_to_keep]
  target_loso_pred_res <- loso_pred_res %>%
    select_("model",
      "region",
      "analysis_time_season",
      "analysis_time_season_week",
      .dots = prediction_cols_to_keep)
  
  return(target_loso_pred_res)
}


#' Fit a stacking model that assigns weights to component models
#' The weights are a function of observed covariates (which?),
#' and are obtained via gradient tree boosting
#' 
#' @param region string with region: either "National" or in the format
#'   "Regionk" for k in {1, ..., 10}
#' @param prediction_target string with either "onset_week", "peak_week",
#'   "peak_inc", "ph1_inc", ..., "ph4_inc"
#' @param component_model_names character vector with names of component models
#' @param loso_preds_path path to directory with leave-one-season-out
#'   predictions from each component model.  Predictions should be saved in
#'   files named like "kde-National-loso-predictions.rds"
#' @param booster what form of boosting to use? see xgboost documentation
#' @param subsample fraction of data to use in bagging.  not supported yet.
#' @param colsample_bytree fraction of explanatory variables to randomly select
#'   in growing each regression tree. see xgboost documentation
#' @param colsample_bylevel fraction of explanatory variables to randomly select
#'   in growing each level of the regression tree. see xgboost documentation
#' @param max_depth maximum depth of regression trees. see xgboost documentation
#' @param min_child_weight not recommended for use. see xgboost documentation
#' @param eta learning rate. see xgboost documentation
#' @param gamma Penalty on number of regression tree leafs. see xgboost documentation
#' @param lambda L2 regularization of contribution to model weights in each
#'   round. see xgboost documentation
#' @param alpha L1 regularization of contribution to model weights in each round.
#'   see xgboost documentation
#' @param nrounds see xgboost documentation
#' @param cv_params optional named list of parameter values to evaluate loss
#'   via cross-validation. Each component is a vector of parameter values with
#'   name one of "booster", "subsample", "colsample_bytree",
#'   "colsample_bylevel", "max_depth", "min_child_weight", "eta", "gamma",
#'    "lambda", "alpha", "nrounds"
#' @param cv_folds list specifying observation groups to use in cross-validation
#'   each list component is a numeric vector of observation indices.
#' @param cv_nfolds integer specifying the number of cross-validation folds to
#'   use.  if cv_folds was provided, cv_nfolds is ignored.  if cv_folds was not
#'   provided, the data will be randomly partitioned into cv_nfolds groups
#' @param cv_refit character describing which of the models specified by the
#'   values in cv_params to refit using the full data set. Either "best",
#'   "ttest", or "none".
#' @param update an object of class xgbstack to update
#' @param nthread how many threads to use. see xgboost documentation
#' @param verbose how much output to generate along the way. 0 for no logging,
#'   1 for some logging
#' 
#' @return a model stacking fit
fit_stacked_model <- function(
  regions,
  prediction_target,
  component_model_names,
  loso_preds_path,
  booster = "gbtree",
  subsample = 1,
  colsample_bytree = 1,
  colsample_bylevel = 1,
  max_depth = 10,
  min_child_weight = -10^10,
  eta = 0.3,
  gamma = 0,
  lambda = 0,
  alpha = 0,
  nrounds = 10,
  cv_params = NULL,
  cv_folds = NULL,
  cv_nfolds = 10L,
  cv_refit = "ttest",
  update = NULL,
  nthread = NULL,
  verbose = 0) {
  require(xgboost)
  require(xgbstack)
  
  ## Assemble training data
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
  
  ## get cross-validation groups if requested
  if(identical(cv_folds, "leave-one-season-out")) {
    cv_folds <- lapply(unique(target_loso_pred_res$analysis_time_season),
      function(season_val) {
        which(target_loso_pred_res$analysis_time_season == season_val)
      })
  }
  
  ## fit stacking model for weights
  fit_formula <- as.formula(paste0(
    paste(component_model_names, collapse = " + "),
    " ~ ",
    "analysis_time_season_week"))
  xgbstack_fit <- xgbstack(fit_formula,
    data = target_loso_pred_res,
    booster = booster,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    colsample_bylevel = colsample_bylevel,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    eta = eta,
    gamma = gamma,
    lambda = lambda,
    alpha = alpha,
    nrounds = nrounds,
    cv_params = cv_params,
    cv_folds = cv_folds,
    cv_nfolds = cv_nfolds,
    cv_refit = cv_refit,
    update = update,
    nthread = nthread,
    verbose = verbose)
  
  return(xgbstack_fit)
}


#' Estimate "optimal" unregularized weights for a linear combination of
#' predictive models.
#' 
#' This function uses either a grid search or the "degenerate EM" algorithm
#' outlined at http://www.cs.cmu.edu/~roni/11761-f16/Presentations/degenerateEM.pdf
#' 
#' @param component_model_log_scores a data frame or matrix of log scores.
#'   Each column gives log scores for a particular predictive model.
#'   Each row corresponds to one observation
#' @param method character, either "em" for degenerate EM or "grid-search" for
#'   grid search
#' @param tol numeric, if method was "em", stopping tolerance
#' @param grid_size number of grid points for each model to evaluate in a grid
#'   search
#' @param return_type character, if method was "grid-search", either "optimal"
#'   to return only the best weights examined or "all" to return the full grid
#' 
#' @return model weights
fit_unregularized_stacked_model <- function(
  component_model_log_scores,
  method = "em",
  tol = .Machine$double.eps,
  grid_size = 101,
  return_type = "optimal") {
  
  component_model_log_scores <- as.matrix(component_model_log_scores)
  
  if(identical(method, "em")) {
    weights <- rep(1 / ncol(component_model_log_scores), ncol(component_model_log_scores))
    log_weights <- log(weights)
    
    prev_score <- -Inf
    curr_score <- mean(logspace_sum_matrix_rows(sweep(component_model_log_scores, 2, log_weights, `+`)))
    
    while((curr_score - prev_score) > tol[1] && (curr_score - prev_score) / abs(curr_score) > tol[1]) {
      log_weight_update_num <- sweep(component_model_log_scores, 2, log_weights, `+`)
      log_weight_update_denom <- logspace_sum_matrix_rows(log_weight_update_num)
      log_weights <-
        logspace_sum_matrix_rows(t(log_weight_update_num - log_weight_update_denom)) -
          log(nrow(component_model_log_scores))
      
      prev_score <- curr_score
      curr_score <- mean(logspace_sum_matrix_rows(sweep(component_model_log_scores, 2, log_weights, `+`)))
    }
    return(exp(log_weights))
  } else if(identical(method, "em-safe")) {
    weights <- rep(1 / ncol(component_model_log_scores), ncol(component_model_log_scores))
    
    prev_score <- -Inf
    curr_score <- mean(log(exp(component_model_log_scores) %*% matrix(weights)))
    
    while((curr_score - prev_score) / abs(curr_score) > tol[1]) {
      weight_update_num <- exp(component_model_log_scores) *
        matrix(rep(weights, each = nrow(component_model_log_scores)), nrow = nrow(component_model_log_scores))
      weight_update_denom <- apply(weight_update_num, 1, sum)
      weights <- apply(weight_update_num / weight_update_denom, 2, mean)
      
      prev_score <- curr_score
      curr_score <- mean(log(exp(component_model_log_scores) %*% matrix(weights)))
    }
    return(weights)
  } else if(identical(method, "grid-search")) {
    weights_grid <- expand.grid(
      lapply(
        seq_len(ncol(component_model_log_scores) - 1),
        function(i) seq(from = 0, to = 1, length = grid_size)
      )
    )
    weights_grid <- weights_grid[apply(weights_grid, 1, sum) <= 1, ]
    weights_grid <- cbind(weights_grid, 1 - apply(weights_grid, 1, sum)) %>%
      `colnames<-`(colnames(component_model_log_scores)) %>%
      mutate(combined_log_score = NA)
    
    for(weights_ind in seq_len(nrow(weights_grid))) {
      weights_grid$combined_log_score[weights_ind] <-
        (
          exp(component_model_log_scores) %*%
            matrix(unlist(weights_grid[weights_ind, seq_len(ncol(component_model_log_scores))]))
          ) %>%
        log() %>%
        mean()
    }
    
    if(identical(return_type, "optimal")) {
      inds <- which(weights_grid$combined_log_score == max(weights_grid$combined_log_score))
    } else {
      inds <- seq_len(nrow(weights_grid))
    }
    return(weights_grid[inds, seq_len(ncol(component_model_log_scores)), drop = FALSE])
  } else {
    stop("invalid method")
  }
}

