library(lubridate)
library(plyr)
library(dplyr)
library(reshape)
library(kcde)
library(doMC)
library(forecast)

### Get command line arguments
args <- commandArgs(trailingOnly=TRUE)

## Some other arguments that we're not using for now, but might want to use later
## I have left in code below that depends on these values.
max_seasonal_lag <- 0
bw_parameterization <- "diagonal"


## national or regional (which region?)
## in order to region as an argument, we inserted a - instead of a space
## between "Region" and region number, so now undo that.
data_set <- gsub("-", " ", args[1])

## Prediction horizon -- integer number of steps ahead to predict
prediction_horizon <- as.integer(args[2])

## Maximum number of non-seasonal lagged observations to explore using for prediction
max_lag <- as.integer(args[3])

## Include terms capturing seasonality?
seasonality <- as.logical(args[4])

## Initialization rule -- "cov" or "scott" for covariance or scott's rule, respectively
bw_init_rule <- "scott"

## Season left out in stacking estimation, we fit separate models
## with each season in the training data left out (leave one season out),
## predictive performance by each model on the held-out season is used
## to estimate the parameters of the stacking ensemble model, which
## assigns weights to each component model
season_left_out <- args[5]

## use data up to but not including first test season to do estimation
first_test_season <- args[6]

## Path to save results in
save_path <- args[7]


### Load data set and set variables describing how the fit is performed
## Load data for nationally reported influenza like illness
data <- read.csv("data-raw/allflu-cleaned.csv", stringsAsFactors = FALSE)

data$time <- as.POSIXct(data$time)

## subset data to be only the region of interest
data <- data[data$region == data_set,]

## Subset data to do estimation using only data up through first test season
first_ind_test_season <- min(which(data$season == first_test_season))
data <- data[seq_len(first_ind_test_season - 1), , drop = FALSE]


prediction_target_var <- "weighted_ili"
data$log_prediction_target <- log(data[, prediction_target_var])
data$seasonally_differenced_log_prediction_target <- 
  ts(c(rep(NA, 52),
       data$log_prediction_target[seq(from = 53, to = nrow(data))] -
         data$log_prediction_target[seq(from = 1, to = nrow(data) - 52)]),
     frequency = 52)

continuous_var_names <- c(
    paste0(c("weighted_ili", "filtered_weighted_ili"), "_horizon", rep(1:52, each=2)),
    paste0(c("weighted_ili", "filtered_weighted_ili"), "_lag", rep(seq(from = 0, to = max_lag + 52 * max_seasonal_lag), each=2))
)
discrete_var_names <- NULL
predictive_vars <- "box_cox_trans_weighted_ili"
time_var <- "time"

target_kernel_fn <- log_pdtmvn_mode_centered_kernel
target_rkernel_fn <- rlog_pdtmvn_mode_centered_kernel
target_initialize_kernel_params_fn <- initialize_params_log_pdtmvn_kernel
target_get_theta_optim_bounds_fn <- get_theta_optim_bounds_log_pdtmvn_kernel
target_vectorize_kernel_params_fn <- vectorize_params_log_pdtmvn_kernel
target_update_theta_from_vectorized_theta_est_fn <- update_theta_from_vectorized_theta_est_log_pdtmvn_kernel
target_discrete_var_range_fns <- NULL

predictive_kernel_fn <- pdtmvn_kernel
predictive_rkernel_fn <- rpdtmvn_kernel
predictive_initialize_kernel_params_fn <- initialize_params_pdtmvn_kernel
predictive_get_theta_optim_bounds_fn <- get_theta_optim_bounds_pdtmvn_kernel
predictive_vectorize_kernel_params_fn <- vectorize_params_pdtmvn_kernel
predictive_update_theta_from_vectorized_theta_est_fn <- update_theta_from_vectorized_theta_est_pdtmvn_kernel
predictive_discrete_var_range_fns <- NULL

variable_selection_method <- "all_included"
crossval_buffer <- ymd("2010-01-01") - ymd("2009-01-01")


### Initial Box Cox transformation of data
## Estimate lambda parameter of box cox transformation using data not including
## season that was left out.
temp_data <- data$weighted_ili
temp_data[data$season == season_left_out] <- NA
lambda <- BoxCox.lambda(temp_data, method = "loglik")
data$box_cox_trans_weighted_ili <- BoxCox(data$weighted_ili, lambda)

### Assemble control parameters for KCDE estimation process


## List describing kernel components -- depends on the values of
## prediction_horizon, filtering, seasonality, and bw_parameterization
kernel_components <- list()

## sample size for initialize_kernel_params_args
if(identical(bw_init_rule, "cov")) {
    init_sample_size <- 1L
} else {
    init_sample_size <- nrow(data)
}

## If requested, periodic kernel component capturing seasonality
if(seasonality) {
    kernel_components <- c(kernel_components,
        list(list(
            vars_and_offsets = data.frame(var_name = "time_index",
                offset_value = 0L,
                offset_type = "lag",
                combined_name = "time_index_lag0",
                stringsAsFactors = FALSE),
            kernel_fn = periodic_kernel,
            theta_fixed = list(period = pi / 365.2425), # 365.2425 is the mean number of days in a year
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_periodic_kernel,
            initialize_kernel_params_args = list(
                sample_size = init_sample_size
            ),
            get_theta_optim_bounds_fn = get_theta_optim_bounds_periodic_kernel,
            get_theta_optim_bounds_args = NULL,
            vectorize_kernel_params_fn = vectorize_params_periodic_kernel,
            vectorize_kernel_params_args = NULL,
            update_theta_from_vectorized_theta_est_fn =
                update_theta_from_vectorized_theta_est_periodic_kernel,
            update_theta_from_vectorized_theta_est_args = NULL
    )))
}

## Kernel components for observed values of incidence
## First step is setup: create list of data frames specifying groups of
## variables and offsets included in each kernel component
lag_values <- NULL

for(seasonal_lag in seq(from = 0, to = max_seasonal_lag)) {
    lag_values <- c(lag_values,
        seq(from = 0, to = max_lag) + 52 * seasonal_lag)
}

if(identical(bw_parameterization, "diagonal")) {
    ## Separate kernel components for each prediction target variable and
    ## predictive variable
    
    vars_and_offsets_groups <- list()
    
    ## Group of variable names and offsets for prediction target
    new_vars_and_offsets_group <- data.frame(
        var_name = prediction_target_var,
        offset_value = prediction_horizon,
        offset_type = "horizon",
        stringsAsFactors = FALSE
    )
    new_vars_and_offsets_group$combined_name <- paste0(
        new_vars_and_offsets_group$var_name,
        "_",
        new_vars_and_offsets_group$offset_type,
        new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- c(vars_and_offsets_groups,
        list(new_vars_and_offsets_group))
    
    ## Groups of variable names and offsets for lagged predictive variables
    
    for(lag_value in lag_values) {
        for(predictive_var in predictive_vars) {
            new_vars_and_offsets_group <- data.frame(
                var_name = predictive_var,
                offset_value = lag_value,
                offset_type = "lag",
                stringsAsFactors = FALSE
            )
            new_vars_and_offsets_group$combined_name <- paste0(
                new_vars_and_offsets_group$var_name,
                "_",
                new_vars_and_offsets_group$offset_type,
                new_vars_and_offsets_group$offset_value
            )
            vars_and_offsets_groups <- c(vars_and_offsets_groups,
                list(new_vars_and_offsets_group))
        }
    }
} else if(identical(bw_parameterization, "full")) {
    ## One kernel component for prediction target variable and all predictive
    ## variables
    
    ## Prediction target variable
    new_vars_and_offsets_group <- data.frame(
        var_name = prediction_target_var,
        offset_value = prediction_horizon,
        offset_type = "horizon",
        stringsAsFactors = FALSE
    )
    
    ## Lagged prediction target == predictive variables
    for(lag_value in lag_values) {
        for(predictive_var in predictive_vars) {
            ## Else, lagged "raw"/unfiltered observed incidence
            new_vars_and_offsets_group <- rbind(
                new_vars_and_offsets_group,
                data.frame(
                    var_name = predictive_var,
                    offset_value = lag_value,
                    offset_type = "lag",
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    
    ## Add combined_name column and put in a list for further processing below
    new_vars_and_offsets_group$combined_name <- paste0(
        new_vars_and_offsets_group$var_name,
        "_",
        new_vars_and_offsets_group$offset_type,
        new_vars_and_offsets_group$offset_value
    )
    vars_and_offsets_groups <- list(new_vars_and_offsets_group)
} else {
    stop("Invalid bandwidth parameterization")
}

## Second step is to actually append the kernel component descriptions to the
## kernel_components list
kernel_components <- c(kernel_components,
    lapply(vars_and_offsets_groups, function(vars_and_offsets) {
        lower_trunc_bds <- rep(-Inf, nrow(vars_and_offsets))
        names(lower_trunc_bds) <- vars_and_offsets$combined_name
        upper_trunc_bds <- rep(Inf, nrow(vars_and_offsets))
        names(upper_trunc_bds) <- vars_and_offsets$combined_name
        
        discrete_var_range_fns <- NULL

        if(prediction_target_var %in% vars_and_offsets$var_name) {
          return(list(
            vars_and_offsets = vars_and_offsets,
            kernel_fn = target_kernel_fn,
            rkernel_fn = target_rkernel_fn,
            theta_fixed = list(
                parameterization = "bw-chol-decomp",
                continuous_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$combined_name %in% continuous_var_names],
                discrete_vars = vars_and_offsets$combined_name[
                    vars_and_offsets$combined_name %in% discrete_var_names],
                discrete_var_range_fns = target_discrete_var_range_fns,
                lower = lower_trunc_bds,
                upper = upper_trunc_bds,
                validate_in_support = FALSE
            ),
            theta_est = list("bw"),
            initialize_kernel_params_fn = target_initialize_kernel_params_fn,
            initialize_kernel_params_args = list(
                sample_size = init_sample_size
            ),
            get_theta_optim_bounds_fn = target_get_theta_optim_bounds_fn,
            get_theta_optim_bounds_args = NULL,
            vectorize_kernel_params_fn = target_vectorize_kernel_params_fn,
            vectorize_kernel_params_args = NULL,
            update_theta_from_vectorized_theta_est_fn = target_update_theta_from_vectorized_theta_est_fn,
            update_theta_from_vectorized_theta_est_args = NULL
        ))
      } else {
        return(list(
          vars_and_offsets = vars_and_offsets,
          kernel_fn = predictive_kernel_fn,
          rkernel_fn = predictive_rkernel_fn,
          theta_fixed = list(
            parameterization = "bw-chol-decomp",
            continuous_vars = vars_and_offsets$combined_name[
              vars_and_offsets$combined_name %in% continuous_var_names],
            discrete_vars = vars_and_offsets$combined_name[
              vars_and_offsets$combined_name %in% discrete_var_names],
            discrete_var_range_fns = predictive_discrete_var_range_fns,
            lower = lower_trunc_bds,
            upper = upper_trunc_bds,
            validate_in_support = FALSE
          ),
          theta_est = list("bw"),
          initialize_kernel_params_fn = predictive_initialize_kernel_params_fn,
          initialize_kernel_params_args = list(
            sample_size = init_sample_size
          ),
          get_theta_optim_bounds_fn = predictive_get_theta_optim_bounds_fn,
          get_theta_optim_bounds_args = NULL,
          vectorize_kernel_params_fn = predictive_vectorize_kernel_params_fn,
          vectorize_kernel_params_args = NULL,
          update_theta_from_vectorized_theta_est_fn = predictive_update_theta_from_vectorized_theta_est_fn,
          update_theta_from_vectorized_theta_est_args = NULL
        ))
      }
    })
)


## Set up filter_control to do no filtering
filter_control <- NULL

## Assemble kernel_components and filter_control created above,
## along with other parameters controling KCDE definition and estimation
kcde_control <- create_kcde_control(X_names = "time_index",
    y_names = prediction_target_var,
    time_name = time_var,
    prediction_horizons = prediction_horizon,
    filter_control = filter_control,
    kernel_components = kernel_components,
    crossval_buffer = crossval_buffer,
    prediction_inds_not_included = which(data$season == season_left_out),
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_args = list(
        prediction_type = "distribution",
        log = TRUE),
    loss_args = NULL,
    par_cores = 4L,
    variable_selection_method = variable_selection_method)



### Do estimation
case_descriptor <- paste0(
    data_set,
    "-prediction_horizon_", prediction_horizon,
    "-max_lag_", max_lag,
    "-seasonality_", seasonality,
    "-season_left_out_", gsub("/", "-", season_left_out)
)

init_theta_vector <- NULL
init_phi_vector <- NULL


registerDoMC(cores = kcde_control$par_cores)

## Get the KCDE fit
fit_time <- system.time({
    kcde_fit <- kcde(data = data,
        kcde_control = kcde_control,
        init_theta_vector = init_theta_vector,
        init_phi_vector = init_phi_vector)
})


### Save results
saveRDS(kcde_fit,
    file = file.path(save_path,
        paste0("kcde_fit-",
            case_descriptor,
            ".rds")
    )
)

saveRDS(lambda,
  file = file.path(save_path,
    paste0("box_cox_lambda-",
      case_descriptor,
      ".rds")
  )
)
