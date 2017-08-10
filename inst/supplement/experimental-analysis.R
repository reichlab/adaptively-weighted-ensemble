library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(awes)
library(mosaic)
library(gridExtra)
library(cowplot)
library(nlme)
library(multcomp)

awes_path <- find.package("awes")

region_season_obs_quantities <- flu_data %>%
  select_(.dots = c("region", "season")) %>%
  distinct() %>%
  filter(season %in% paste0(2011:2015, "/", 2012:2016)) %>%
  mutate(
    observed_onset_week = NA,
    observed_peak_week = NA
  )

for(rs_row in seq_len(nrow(region_season_obs_quantities))) {
  temp <- get_observed_seasonal_quantities(
    data = flu_data[flu_data$region == region_season_obs_quantities$region[rs_row], , drop = FALSE],
    season = region_season_obs_quantities$season[rs_row],
    first_CDC_season_week = 10,
    last_CDC_season_week = 42,
    onset_baseline =
      get_onset_baseline(region = region_season_obs_quantities$region[rs_row],
        season = region_season_obs_quantities$season[rs_row]),
    incidence_var = "weighted_ili",
    incidence_bins = data.frame(
      lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
      upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
    incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
  )
  
  region_season_obs_quantities$observed_onset_week[rs_row] <-
    temp$observed_onset_week
  region_season_obs_quantities$observed_peak_week[rs_row] <-
    temp$observed_peak_week[1]
}

region_season_obs_quantities$observed_onset_week[
  region_season_obs_quantities$observed_onset_week == "none"] <- 42
region_season_obs_quantities <- region_season_obs_quantities %>%
  transmute(
    region = ifelse(region == "X", "National", gsub(" ", "", region)),
    analysis_time_season = season,
    observed_onset_week = as.numeric(observed_onset_week),
    observed_peak_week = observed_peak_week
  )

all_models <- c("kde", "kcde", "sarima", "equal_weights", "em_stacking", "xgb_stacking_unregularized", "xgb_stacking_reg_w", "xgb_stacking_reg_wu", "xgb_stacking_reg_wui")
all_targets <- c("onset", "peak_week", "peak_inc")
preds <- assemble_predictions(
  preds_path = file.path(awes_path, "evaluation/test-predictions"),
  models = all_models,
  prediction_targets = all_targets,
  prediction_types = "log_score"
) %>%
  filter(analysis_time_season_week %in% 10:40) %>%
  gather_("prediction_target", "log_score",
    c("onset_log_score", "peak_week_log_score", "peak_inc_log_score")) %>%
  mutate(
    prediction_target = substr(prediction_target, 1, nchar(prediction_target) - 10),
    score = exp(log_score)
  )

preds$log_score[is.infinite(preds$log_score)] <- -10
preds <- preds %>%
  left_join(region_season_obs_quantities, by = c("region", "analysis_time_season")) %>%
  mutate(
    before_onset = (analysis_time_season_week < observed_onset_week),
    before_peak = (analysis_time_season_week < observed_peak_week))

preds$before_target_date <- ifelse(
  preds$prediction_target == "onset",
  c("on or after\ntarget date", "before\ntarget date")[preds$before_onset + 1],
  c("on or after\ntarget date", "before\ntarget date")[preds$before_peak + 1]
)

preds$model[preds$model == "kde"] <- "KDE"
preds$model[preds$model == "kcde"] <- "KCDE"
preds$model[preds$model == "sarima"] <- "SARIMA"
preds$model[preds$model == "equal_weights"] <- "EW"
preds$model[preds$model == "em_stacking"] <- "CW"
preds$model[preds$model == "xgb_stacking_unregularized"] <- "FW-wu"
preds$model[preds$model == "xgb_stacking_reg_w"] <- "FW-reg-w"
preds$model[preds$model == "xgb_stacking_reg_wu"] <- "FW-reg-wu"
preds$model[preds$model == "xgb_stacking_reg_wui"] <- "FW-reg-wui"
preds$model <- factor(preds$model,
  levels = c("KDE", "KCDE", "SARIMA",
    "EW", "CW", "FW-wu",
    "FW-reg-w", "FW-reg-wu", "FW-reg-wui"))

preds$prediction_target[preds$prediction_target == "onset"] <- "Onset Timing"
preds$prediction_target[preds$prediction_target == "peak_inc"] <- "Peak Incidence"
preds$prediction_target[preds$prediction_target == "peak_week"] <- "Peak Timing"
preds$prediction_target <- factor(preds$prediction_target,
  levels = c("Onset Timing", "Peak Timing", "Peak Incidence"))

preds <- preds %>%
  mutate(
    region = factor(region),
    analysis_time_season = factor(analysis_time_season)
  )

preds_median_ls_by_model_target <- preds %>%
  group_by_("model", "prediction_target") %>%
  summarize(ls_median = median(log_score))

preds <- preds %>%
  left_join(preds_median_ls_by_model_target) %>%
  mutate(abs_log_score_diff_from_median = abs(log_score - ls_median))

ggplot(data = preds) +
  geom_density(aes(x = abs_log_score_diff_from_median, color = prediction_target)) +
  facet_wrap( ~ model) +
  theme_bw()

preds %>%
  group_by(model, prediction_target) %>%
  summarize(mean_diff = mean(abs_log_score_diff_from_median)) %>%
  arrange(prediction_target, mean_diff) %>%
  as.data.frame()

aov_fit <- lm(abs_log_score_diff_from_median ~ model,
  data = preds %>% filter(prediction_target == "Peak Timing"))
anova(aov_fit)
summary(aov_fit)

unique_region_season_model_autocor <- preds %>%
#  filter(before_target_date == "before\ntarget date") %>%
  select_("region", "analysis_time_season", "model", "prediction_target", "abs_log_score_diff_from_median") %>%
  group_by_("region", "analysis_time_season", "model", "prediction_target") %>%
  summarize(
    log_score_autocor = acf(abs_log_score_diff_from_median, plot = FALSE)[["acf"]][2, 1, 1]
  )
## set to 1 for kde
unique_region_season_model_autocor$log_score_autocor[
  is.na(unique_region_season_model_autocor$log_score_autocor)] <- 1

preds <- preds %>%
  mutate(
    region = factor(region),
    analysis_time_season = factor(analysis_time_season),
    model = factor(model),
    prediction_target = factor(prediction_target),
    unique_region_season_model_target = factor(
      paste(region, analysis_time_season, model, prediction_target, sep = "_")
    ),
    unique_model_target = factor(paste(model, prediction_target, sep = "_")),
    unique_region_season_model = factor(
      paste(region, analysis_time_season, model, sep = "_")
    )
  )

reduced_lm_fit <- lm(mean_diff ~ model * prediction_target,
  data = preds %>%
    mutate(model = factor(as.character(model),
      levels = c("FW-reg-w", "KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-wu", "FW-reg-wui"))) %>%
    group_by(model, prediction_target, analysis_time_season, region) %>%
    summarize(mean_diff = mean(abs_log_score_diff_from_median)))

preds %>%
  group_by(model, prediction_target, analysis_time_season, region) %>%
  summarize(mean_diff = mean(abs_log_score_diff_from_median))

preds %>%
  group_by(model, prediction_target) %>%
  summarize(mean_diff = mean(abs_log_score_diff_from_median)) %>%
  arrange(prediction_target, mean_diff)


gls_fit <- gls(abs_log_score_diff_from_median ~ model * prediction_target,
  data=preds,
  correlation=corARMA(p=1, q=0, form=~analysis_time_season_week | unique_region_season_model_target))

gls_fit_onset <- gls(  abs_log_score_diff_from_median ~ model,
  data = preds %>%
    filter(prediction_target == "Onset Timing") %>%
    mutate(model = factor(as.character(model),
      levels = c("EW", "KDE", "KCDE", "SARIMA", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui"))),
  correlation=corARMA(p=1, q=0, form=~analysis_time_season_week | unique_region_season_model))

gls_fit_pt <- gls(  abs_log_score_diff_from_median ~ model,
  data = preds %>%
    filter(prediction_target == "Peak Timing") %>%
    mutate(model = factor(as.character(model),
      levels = c("EW", "KDE", "KCDE", "SARIMA", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui"))),
  correlation=corARMA(p=1, q=0, form=~analysis_time_season_week | unique_region_season_model))

gls_fit_pi <- gls(  abs_log_score_diff_from_median ~ model,
  data = preds %>%
    filter(prediction_target == "Peak Incidence") %>%
    mutate(model = factor(as.character(model),
      levels = c("EW", "KDE", "KCDE", "SARIMA", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui"))),
  correlation=corARMA(p=1, q=0, form=~analysis_time_season_week | unique_region_season_model_target))


lm_fit <- lme(
  abs_log_score_diff_from_median ~ model * prediction_target,
#  random = ~ 1 | unique_region_season_model_target,
  random = ~ 1 | unique_model_target,
  correlation = corAR1(
    value = mean(unique_region_season_model_autocor$log_score_autocor)#,
#    form = ~ analysis_time_season_week | unique_region_season_model_target),
  ),
  method = "REML",
  data = preds)
#  data = preds[preds$before_target_date == "before\ntarget date", ])

all_models <- as.character(unique(preds$model))
all_targets <- as.character(unique(preds$prediction_target))
num_models <- length(all_models)
num_targets <- length(all_targets)

unique_model_descriptors <- paste0("model", all_models)
unique_target_descriptors <- paste0("prediction_target", all_targets)

lc_df <- expand.grid(
  model = all_models,
  target = all_targets,
  stringsAsFactors = FALSE)

lc_df$model_descriptor <- paste0("model", lc_df$model)
lc_df$target_descriptor <- paste0("prediction_target", lc_df$target)
lc_df$name <- apply(as.matrix(lc_df[, 1:2]), 1, paste, collapse = "-")

num_leading_cols <- ncol(lc_df)
coef_cols <- seq(
  from = num_leading_cols + 1,
  length = num_models * num_targets
)

# corresponding indicator vector for each coefficient
#coef_names <- names(fixef(lm_fit))
coef_names <- names(coef(gls_fit))
unique_coef_name_component_descriptors <- unique(unlist(strsplit(coef_names, ":")))
intercept_model <- unique_model_descriptors[
  !(unique_model_descriptors %in% unique_coef_name_component_descriptors)]
intercept_target <- unique_target_descriptors[
  !(unique_target_descriptors %in% unique_coef_name_component_descriptors)]
for(coef_ind in seq(from = 1, to = length(coef_names))) {
  split_name <- unlist(strsplit(coef_names[[coef_ind]], ":"))
  if(!any(split_name %in% unique_model_descriptors[unique_model_descriptors != intercept_model])) {
    split_name <- c(split_name, unique_model_descriptors)
  }
  if(!any(split_name %in% unique_target_descriptors[unique_target_descriptors != intercept_target])) {
    split_name <- c(split_name, unique_target_descriptors)
  }
  
  lc_df[[paste0("coef", coef_ind)]] <- 0
  lc_df[[paste0("coef", coef_ind)]][
    lc_df$model_descriptor %in% split_name &
      lc_df$target_descriptor %in% split_name] <- 1
}

model_row_inds <- seq_len(nrow(lc_df))

# # ## contrasts of (mean performance model 1) - (mean performance model 2) for all model pairs
# rowind <- nrow(lc_df)
# for(prediction_target in all_targets) {
#   for(fit_method1_ind in seq_len(length(all_models) - 1)) {
#     for(fit_method2_ind in seq(from = fit_method1_ind + 1, to = length(all_models))) {
#       rowind <- rowind + 1
#       
#       m1_name <- all_models[fit_method1_ind]
#       m2_name <- all_models[fit_method2_ind]
#       
#       m1_rowind <- which(lc_df$name == paste(m1_name, prediction_target, sep = "-"))
#       m2_rowind <- which(lc_df$name == paste(m2_name, prediction_target, sep = "-"))
#       
#       lc_df[rowind, ] <- rep(NA, ncol(lc_df))
#       lc_df$name[rowind] <- paste0(m1_name, "-", m2_name, "-", prediction_target)
#       lc_df$target[rowind] <- prediction_target
#       lc_df[rowind, coef_cols] <- lc_df[m1_rowind, coef_cols] - lc_df[m2_rowind, coef_cols]
#     }
#   }
# }
# 
# contrast_row_inds <- seq(from = length(model_row_inds) + 1, to = nrow(lc_df))

## average performance of each model across all three prediction targets
for(fit_method_ind in seq_len(length(all_models))) {
  rowind <- rowind + 1
  
  m_name <- all_models[fit_method_ind]
  
  m_rowinds <- which(lc_df$model == m_name)
  
  lc_df[rowind, ] <- rep(NA, ncol(lc_df))
  lc_df$model[rowind] <- m_name
  lc_df$name[rowind] <- paste0(m_name, "-grand_mean")
  lc_df[rowind, coef_cols] <- apply(lc_df[m_rowinds, coef_cols], 2, mean)
}

#lc_df <- lc_df[c(1:9, 136:144), ]

lc_df$name <- factor(lc_df$name, levels = unique(lc_df$name))

K_mat <- as.matrix(lc_df[, coef_cols])

# get point estimates
#lc_df$pt_est <- as.vector(K_mat %*% matrix(fixef(lm_fit)))
lc_df$pt_est <- as.vector(K_mat %*% matrix(coef(gls_fit)))

# get familywise CIs
confint_rows <- seq_len(nrow(lc_df))
lc_df$fam_CI_lb <- NA
lc_df$fam_CI_ub <- NA
#fam_CI_obj <- glht(lm_fit, linfct = K_mat[confint_rows, ])
fam_CI_obj <- glht(gls_fit, linfct = K_mat[confint_rows, ])
temp <- confint(fam_CI_obj)$confint
lc_df$fam_CI_lb[confint_rows] <- temp[, 2]
lc_df$fam_CI_ub[confint_rows] <- temp[, 3]

# get individual CIs
lc_df$ind_CI_lb <- NA
lc_df$ind_CI_ub <- NA
for(rowind in confint_rows) {
#  ind_CI_obj <- glht(lm_fit, linfct = K_mat[rowind, , drop = FALSE])
  ind_CI_obj <- glht(gls_fit, linfct = K_mat[rowind, , drop = FALSE])
  temp <- confint(ind_CI_obj)$confint
  lc_df$ind_CI_lb[rowind] <- temp[, 2]
  lc_df$ind_CI_ub[rowind] <- temp[, 3]
}























model_pair <- c("FW-reg-w", "SARIMA")
gls_fit <- gls(abs_log_score_diff_from_median ~ model * prediction_target,
  data = preds %>%
    filter(model %in% model_pair),
  correlation =
    corARMA(p=1, q=0, form=~analysis_time_season_week | unique_region_season_model_target))

all_models <- model_pair
all_targets <- as.character(unique(preds$prediction_target))
num_models <- length(all_models)
num_targets <- length(all_targets)

unique_model_descriptors <- paste0("model", all_models)
unique_target_descriptors <- paste0("prediction_target", all_targets)

lc_df <- expand.grid(
  model = all_models,
  target = all_targets,
  stringsAsFactors = FALSE)

lc_df$model_descriptor <- paste0("model", lc_df$model)
lc_df$target_descriptor <- paste0("prediction_target", lc_df$target)
lc_df$name <- apply(as.matrix(lc_df[, 1:2]), 1, paste, collapse = "-")

num_leading_cols <- ncol(lc_df)
coef_cols <- seq(
  from = num_leading_cols + 1,
  length = num_models * num_targets
)

# corresponding indicator vector for each coefficient
#coef_names <- names(fixef(lm_fit))
coef_names <- names(coef(gls_fit))
unique_coef_name_component_descriptors <- unique(unlist(strsplit(coef_names, ":")))
intercept_model <- unique_model_descriptors[
  !(unique_model_descriptors %in% unique_coef_name_component_descriptors)]
intercept_target <- unique_target_descriptors[
  !(unique_target_descriptors %in% unique_coef_name_component_descriptors)]
for(coef_ind in seq(from = 1, to = length(coef_names))) {
  split_name <- unlist(strsplit(coef_names[[coef_ind]], ":"))
  if(!any(split_name %in% unique_model_descriptors[unique_model_descriptors != intercept_model])) {
    split_name <- c(split_name, unique_model_descriptors)
  }
  if(!any(split_name %in% unique_target_descriptors[unique_target_descriptors != intercept_target])) {
    split_name <- c(split_name, unique_target_descriptors)
  }
  
  lc_df[[paste0("coef", coef_ind)]] <- 0
  lc_df[[paste0("coef", coef_ind)]][
    lc_df$model_descriptor %in% split_name &
      lc_df$target_descriptor %in% split_name] <- 1
}

model_row_inds <- seq_len(nrow(lc_df))

# # ## contrasts of (mean performance model 1) - (mean performance model 2) for all model pairs
rowind <- nrow(lc_df)
# for(prediction_target in all_targets) {
#   for(fit_method1_ind in seq_len(length(all_models) - 1)) {
#     for(fit_method2_ind in seq(from = fit_method1_ind + 1, to = length(all_models))) {
#       rowind <- rowind + 1
#       
#       m1_name <- all_models[fit_method1_ind]
#       m2_name <- all_models[fit_method2_ind]
#       
#       m1_rowind <- which(lc_df$name == paste(m1_name, prediction_target, sep = "-"))
#       m2_rowind <- which(lc_df$name == paste(m2_name, prediction_target, sep = "-"))
#       
#       lc_df[rowind, ] <- rep(NA, ncol(lc_df))
#       lc_df$name[rowind] <- paste0(m1_name, "-", m2_name, "-", prediction_target)
#       lc_df$target[rowind] <- prediction_target
#       lc_df[rowind, coef_cols] <- lc_df[m1_rowind, coef_cols] - lc_df[m2_rowind, coef_cols]
#     }
#   }
# }
# 
# contrast_row_inds <- seq(from = length(model_row_inds) + 1, to = nrow(lc_df))

## average performance of each model across all three prediction targets
for(fit_method_ind in seq_len(length(all_models))) {
  rowind <- rowind + 1
  
  m_name <- all_models[fit_method_ind]
  
  m_rowinds <- which(lc_df$model == m_name)
  
  lc_df[rowind, ] <- rep(NA, ncol(lc_df))
  lc_df$model[rowind] <- m_name
  lc_df$name[rowind] <- paste0(m_name, "-grand_mean")
  lc_df[rowind, coef_cols] <- apply(lc_df[m_rowinds, coef_cols], 2, mean)
}

## difference in average performance across all three prediction targets between the two methods
rowind <- rowind + 1
m1_name <- all_models[1]
m2_name <- all_models[2]
m1_rowind <- which(lc_df$name == paste0(m1_name, "-grand_mean"))
m2_rowind <- which(lc_df$name == paste0(m2_name, "-grand_mean"))
lc_df[rowind, ] <- rep(NA, ncol(lc_df))
lc_df$name[rowind] <- paste0(m1_name, "-", m2_name, "-mean")
lc_df[rowind, coef_cols] <- lc_df[m1_rowind, coef_cols] - lc_df[m2_rowind, coef_cols]

lc_df <- lc_df[rowind, , drop = FALSE]

#lc_df <- lc_df[c(1:9, 136:144), ]

lc_df$name <- factor(lc_df$name, levels = unique(lc_df$name))

K_mat <- as.matrix(lc_df[, coef_cols])

# get point estimates
#lc_df$pt_est <- as.vector(K_mat %*% matrix(fixef(lm_fit)))
lc_df$pt_est <- as.vector(K_mat %*% matrix(coef(gls_fit)))

# get familywise CIs
confint_rows <- seq_len(nrow(lc_df))
lc_df$fam_CI_lb <- NA
lc_df$fam_CI_ub <- NA
#fam_CI_obj <- glht(lm_fit, linfct = K_mat[confint_rows, ])
fam_CI_obj <- glht(gls_fit, linfct = K_mat[confint_rows, , drop = FALSE])
temp <- confint(fam_CI_obj)$confint
lc_df$fam_CI_lb[confint_rows] <- temp[, 2]
lc_df$fam_CI_ub[confint_rows] <- temp[, 3]

# get individual CIs
lc_df$ind_CI_lb <- NA
lc_df$ind_CI_ub <- NA
for(rowind in confint_rows) {
  #  ind_CI_obj <- glht(lm_fit, linfct = K_mat[rowind, , drop = FALSE])
  ind_CI_obj <- glht(gls_fit, linfct = K_mat[rowind, , drop = FALSE])
  temp <- confint(ind_CI_obj)$confint
  lc_df$ind_CI_lb[rowind] <- temp[, 2]
  lc_df$ind_CI_ub[rowind] <- temp[, 3]
}








preds_wide <- preds %>%
  dplyr::filter(before_target_date == "before\ntarget date") %>%
  dplyr::select(-score) %>%
  tidyr::spread(model, log_score)


summarized_res <-
  preds_wide[, c("KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui", "prediction_target")] %>%
  #  select_(.dots = c("KDE", "KCDE", "SARIMA", "EW", "CW", `FW-wu`, `FW-reg-w`, `FW-reg-wu`, `FW-reg-wui`, "prediction_target")) %>%
  #  select_(.dots = c("KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui", "prediction_target")) %>%
  group_by(prediction_target) %>%
  summarize_each(funs(mean, var)) %>%
  t() %>%
  `colnames<-`(.[1, ])
for(model_val in c("KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui")) {
  cat(model_val)
  for(target_val in colnames(summarized_res)) {
    cat(" & ")
    cat(format(round(as.numeric(summarized_res[paste0(model_val, "_mean"), target_val]), 2), nsmall = 2))
    cat(" & ")
    cat(format(round(as.numeric(summarized_res[paste0(model_val, "_var"), target_val]), 2), nsmall = 2))
    if(target_val %in% colnames(summarized_res)[1:2]) {
      cat(" & ")
    }
  }
  cat(" \\\\\n")
}
