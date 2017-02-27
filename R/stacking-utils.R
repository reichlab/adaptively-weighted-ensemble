## utility functions that might be useful for multiple prediction methods

#' Compute and save a weighted combination of prediction results for the test
#' phase.
#' 
#' @param preds_path path to test phase prediction results
#' @param region character specifying a single region in the form of "National"
#'   or "Regionk" for k in 1, ..., 10
#' @param models character vector specifying component models whose predictions
#'   should be combined
#' @param matrix of model weights.  one column for each entry in "models", in
#'   the same order.  one row for each row in 
#' @param weighting_method character string specifying method name for
#'   weighting; used as part of save file name
#' 
#' @return weighted combination of predictions
#' 
#' @export
weighted_combine_predictions <- function(
  component_preds,
  weights) {
  if(length(region) != 1L) {
    stop("region must be length 1")
  }
  component_preds <- assemble_predictions(
    preds_path = preds_path,
    regions = region,
    models = models,
    prediction_targets = prediction_target
  )
  
  
}


#' Assemble leave-one-season-out or test phase predictions made by kde, kcde, and sarima
#' models on training data.
#' 
#' @param preds_path path to leave-one-season-out or test phase predictions.
#' @param regions character vector specifying regions for which to get predictions,
#'   "National" or "Regionk" for k in 1, ..., 10.
#' @param seasons character vector specifying seasons for which to get predictions,
#'   "2011/2012" or "2011-2012" or "*" for all seasons
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
assemble_predictions <- function(
  preds_path = "inst/estimation/loso-predictions",
  regions = c("National", paste0("Region", 1:10)),
  seasons = "*",
  models = c("kde", "kcde", "sarima"),
  prediction_targets = c("onset", "peak_week", "peak_inc", "ph_1_inc", "ph_2_inc", "ph_3_inc", "ph_4_inc"),
  prediction_types = c("log_score", "competition_log_score", "bin_log_probs")
  ) {
  ## Seasons in format used in prediction file names
  seasons <- gsub("/", "-", seasons)
  
  ## names of files with prediction results to load
  model_region_season_combos <- expand.grid(
    model = models,
    region = regions,
    season = seasons,
    stringsAsFactors = TRUE
  )
  file_name_patterns <- apply(model_region_season_combos, 1,
    function(combo) {
      paste0(preds_path, "/", combo["model"], "-", combo["region"], "-", combo["season"], "*")
    })
  file_names <- Sys.glob(file_name_patterns)
  
  ## load the prediction results
  pred_res <- rbind.fill(lapply(
    file_names,
    function(file_path) {
      region_val <- names(unlist(sapply(regions, function(region_val) grep(paste0(region_val, "-"), file_path))))
      readRDS(file_path) %>%
        mutate(region = region_val)
    }
  ))
  
  ## narrow down to the specified prediction targets and types
  prediction_cols_to_keep_templates <- outer(prediction_targets, prediction_types,
    function(target, type) {
      paste0(target, "_", ifelse(type == "bin_log_probs", "bin_.*_log_prob", type))
    }) %>%
    as.vector()
  prediction_cols_to_keep <- lapply(
    prediction_cols_to_keep_templates,
    function(pattern) grep(pattern, names(pred_res))) %>%
    unlist()
  prediction_cols_to_keep <- names(pred_res)[prediction_cols_to_keep]
  target_pred_res <- pred_res %>%
    select_("model",
      "region",
      "analysis_time_season",
      "analysis_time_season_week",
      .dots = prediction_cols_to_keep)
  
  return(target_pred_res)
}

#' Assemble a data frame of inputs to stacking
#' 
#' @param regions string with region: either "National" or in the format
#'   "Regionk" for k in {1, ..., 10}
#' @param prediction_target string with either "onset", "peak_week",
#'   "peak_inc", "ph1_inc", ..., "ph4_inc"
#' @param component_model_names character vector with names of component models
#' @param explanatory_variables character vector with names of explanatory variables
#'   to include for weights; a non-empty subset of "analysis_time_season_week",
#'   "kcde_model_confidence", "sarima_model_confidence", "weighted_ili"
#' @param include_model_performance boolean; should measures of model
#'   performance be included in the return result? Generally, should be TRUE if
#'   we're getting inputs for model training and FALSE if we're just getting
#'   covariates to calculate model weights.
#' @param preds_path path to directory with leave-one-season-out or test phase
#'   predictions from each component model.  Predictions should be saved in
#'   files named like "kde-National-loso-predictions.rds"
#'
#' @return a data frame with covariates that model weights depend on
#'   (as specified by explanatory_variables) and possibly log scores from each
#'   component model
#' 
#' @export
assemble_stacking_inputs <- function(
  regions,
  prediction_target,
  component_model_names,
  explanatory_variables,
  include_model_performance = FALSE,
  preds_path
) {
  ## Load prediction results
  target_pred_res <- assemble_predictions(
    preds_path = preds_path,
    regions = regions,
    models = component_model_names,
    prediction_targets = prediction_target
  )
  
  ## get model confidence
  cols_to_examine <- grep(paste0(prediction_target, ".*_log_prob"), colnames(target_pred_res), value = TRUE)
  target_pred_res$model_confidence <- sapply(
    seq_len(nrow(target_pred_res)),
    function(i) {
      temp <- cumsum(sort(unlist(exp(target_pred_res[i, cols_to_examine])), decreasing = TRUE))
      temp <- temp / temp[length(temp)]
      min(which(temp >= 0.90))
    })
  target_pred_res_with_model_confidence <- target_pred_res %>%
    select_(.dots = c("region", "analysis_time_season", "analysis_time_season_week", "model", "model_confidence")) %>%
    spread_("model", "model_confidence")
  colnames(target_pred_res_with_model_confidence) <- c("region",
    "analysis_time_season",
    "analysis_time_season_week",
    paste0(colnames(target_pred_res_with_model_confidence)[
      seq(from = ncol(target_pred_res_with_model_confidence) - 2, to = ncol(target_pred_res_with_model_confidence))], "_model_confidence"))
  
  if(include_model_performance) {
    ## spread log scores
    target_pred_res_with_log_score <- target_pred_res %>%
      select_(.dots = c("region", "analysis_time_season", "analysis_time_season_week", "model", paste0(prediction_target, "_log_score"))) %>%
      spread_("model", paste0(prediction_target, "_log_score"))
    colnames(target_pred_res_with_log_score) <- c("region",
      "analysis_time_season",
      "analysis_time_season_week",
      paste0(colnames(target_pred_res_with_log_score)[
        seq(from = ncol(target_pred_res_with_log_score) - 2, to = ncol(target_pred_res_with_log_score))], "_log_score"))
    
    ## join to target_pred_res_with_log_score and store in target_pred_res
    target_pred_res <- left_join(target_pred_res_with_log_score,
      target_pred_res_with_model_confidence,
      by = c("region",
        "analysis_time_season",
        "analysis_time_season_week")) %>%
      mutate(region = ifelse(region == "National", "National", paste0("Region ", substr(region, 7, nchar(region)))))
  } else {
    target_pred_res <- target_pred_res_with_model_confidence %>%
      mutate(region = ifelse(region == "National", "National", paste0("Region ", substr(region, 7, nchar(region)))))
  }
  
  ## add weighted ili
  ## Load and clean up data set

  ## subset data to be only the region of interest
  regions_data <- sapply(regions,
    function(reg) {
      if(reg == "National") {
        return("X")
      } else {
        return(paste0("Region ", substr(reg, 7, nchar(reg))))
      }
    })
  data <- flu_data[flu_data$region %in% regions_data,]
  
  ## join to target_pred_res
  reduced_data <- data %>%
    as.data.frame() %>%
    transmute(
      region = ifelse(region == "X", "National", as.character(region)),
      analysis_time_season = season,
      analysis_time_season_week = season_week,
      weighted_ili = weighted_ili
    )
  target_pred_res <- left_join(target_pred_res,
    reduced_data,
    by = c("region",
      "analysis_time_season",
      "analysis_time_season_week"))
  
  ## drop rows where any of the (response or) explanatory variables for weights
  ## have NA values.  this is aggressive.
  if(include_model_performance) {
    vars_to_check <- c(paste0(component_model_names, "_log_score"), explanatory_variables)
  } else {
    vars_to_check <- explanatory_variables
  }
  target_pred_res <- target_pred_res[
    apply(target_pred_res[, vars_to_check, drop = FALSE],
      1,
      function(x) {!any(is.na(x))}), # row has na?  drop if so
    , # all columns
    drop = FALSE]
  
  return(target_pred_res)
}


#' Fit a stacking model that assigns weights to component models
#' The weights are a function of observed covariates (which?),
#' and are obtained via gradient tree boosting
#' 
#' @param regions string with region: either "National" or in the format
#'   "Regionk" for k in {1, ..., 10}
#' @param prediction_target string with either "onset", "peak_week",
#'   "peak_inc", "ph1_inc", ..., "ph4_inc"
#' @param component_model_names character vector with names of component models
#' @param explanatory_variables character vector with names of explanatory variables
#'   to include for weights; a non-empty subset of "analysis_time_season_week",
#'   "kcde_model_confidence", "sarima_model_confidence", "weighted_ili"
#' @param loso_preds_path path to directory with leave-one-season-out
#'   predictions from each component model.  Predictions should be saved in
#'   files named like "kde-National-loso-predictions.rds"
#' @param seasons_to_leave_out optional character vector of seasons to leave out
#'   of stacking estimation
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
#' 
#' @export
fit_stacked_model <- function(
  regions,
  prediction_target,
  component_model_names,
  explanatory_variables =
    c("analysis_time_season_week", "kcde_model_confidence", "sarima_model_confidence", "weighted_ili"),
  loso_preds_path,
  seasons_to_leave_out,
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
  
  ### Assemble training data
  target_loso_pred_res <- assemble_stacking_inputs(
    regions = regions,
    prediction_target = prediction_target,
    component_model_names = component_model_names,
    explanatory_variables = explanatory_variables,
    include_model_performance = TRUE,
    preds_path = loso_preds_path
  )
  
  ## if requested, drop season(s)
  if(!missing(seasons_to_leave_out)) {
    target_loso_pred_res <- 
      target_loso_pred_res[!(target_loso_pred_res$analysis_time_season %in% seasons_to_leave_out), ]
  }
  
  ## get cross-validation groups if requested
  if(identical(cv_folds, "leave-one-season-out")) {
    cv_folds <- lapply(unique(target_loso_pred_res$analysis_time_season),
      function(season_val) {
        which(target_loso_pred_res$analysis_time_season == season_val)
      })
  }
  
  ## fit stacking model for weights
  fit_formula <- as.formula(paste0(
    paste0(component_model_names, "_log_score", collapse = " + "),
    " ~ ",
    paste(explanatory_variables, collapse = " + ")))
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
#' 
#' @export
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


#' Estimate weights for a linear combination of predictive models.
#' 
#' This function uses the "degenerate EM" algorithm outlined at
#' http://www.cs.cmu.edu/~roni/11761-f16/Presentations/degenerateEM.pdf,
#' but modified in a way that is ad hoc but I think could be justified
#' to use kernel-weighted observations.
#' 
#' @param component_model_log_scores a data frame or matrix of log scores.
#'   Each column gives log scores for a particular predictive model.
#'   Each row corresponds to one observation
#' @param covariate a single covariate over which weights should be smoothed
#'   (more than one covariate may be supported later)
#' @param bw a bandwidth for the smoothing
#' @param tol numeric, if method was "em", stopping tolerance
#' 
#' @return model weights
#' 
#' @export
fit_kernel_smoothed_stacked_model_fixed_bw <- function(
  component_model_log_scores,
  covariate,
  prediction_covariate,
  bw,
  tol = .Machine$double.eps) {
  
  component_model_log_scores <- as.matrix(component_model_log_scores)
  covariate <- as.matrix(covariate)
  prediction_covariate <- as.matrix(prediction_covariate)
  method <- "em-safe"
  
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
    combined_weights <- matrix(NA, nrow = nrow(prediction_covariate),
      ncol = ncol(component_model_log_scores))
    for(i in seq_len(nrow(combined_weights))) {
      cat(paste0("Getting weights for covariate ", i, " of ", nrow(combined_weights), ".\n"))
      cat(Sys.time())
      cat("\n")
      obs_weights <- dnorm(covariate, prediction_covariate[i, ], sd = bw)
      obs_weights <- obs_weights / sum(obs_weights)
      
      weights <- rep(1 / ncol(component_model_log_scores), ncol(component_model_log_scores))
      
      prev_score <- -Inf
      curr_score <- weighted.mean(
        log(exp(component_model_log_scores) %*% matrix(weights)),
        w = obs_weights)
      
      while((curr_score - prev_score) / abs(curr_score) > tol[1]) {
        weight_update_num <- exp(component_model_log_scores) *
          matrix(rep(weights, each = nrow(component_model_log_scores)), nrow = nrow(component_model_log_scores))
        weight_update_denom <- apply(weight_update_num, 1, sum)
        weights <- apply(weight_update_num / weight_update_denom, 2, weighted.mean, w = obs_weights)
        
        prev_score <- curr_score
        curr_score <- weighted.mean(
          log(exp(component_model_log_scores) %*% matrix(weights)),
          w = obs_weights)
#        cat(curr_score)
      }
      combined_weights[i, ] <- weights
    }
    return(combined_weights)
  } else {
    stop("invalid method")
  }
}

