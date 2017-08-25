library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(awes)
library(nlme)
library(multcomp)
library(doMC)
library(mosaic)

awes_path <- find.package("awes")

all_models <- c("KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui")

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

preds$log_score[is.infinite(preds$log_score)] <- -15
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


preds_wide <- preds %>%
#  dplyr::filter(before_target_date == "before\ntarget date") %>%
  dplyr::select(-score) %>%
  tidyr::spread(model, log_score)

registerDoMC(cores = 4)

all_models <- c("KDE", "KCDE", "SARIMA", "EW", "CW", "FW-wu", "FW-reg-w", "FW-reg-wu", "FW-reg-wui")
m_combos <- expand.grid(
  m1_ind = seq_along(all_models),
  m2_ind = seq_along(all_models)) %>%
  filter(m1_ind < m2_ind)
m_combos$m1 <- all_models[m_combos$m1_ind]
m_combos$m2 <- all_models[m_combos$m2_ind]

get_one_perm <- function(preds_wide_perm, all_models, group_by_vars = NULL) {
  model_cols <- which(colnames(preds_wide_perm) %in% all_models)

  model_scores <- preds_wide_perm[, model_cols]
  model_cols <- seq_along(model_cols)

  if(is.null(group_by_vars)) {
    for(i in seq_len(nrow(preds_wide_perm))) {
      model_scores[i, model_cols] <- base::sample(model_scores[i, model_cols])
    }
  } else {
    groups <- preds_wide_perm[, group_by_vars] %>%
      distinct() %>%
      as.data.frame()
    for(group_ind in seq_len(nrow(groups))) {
      preds_inds <- which(
        preds_wide_perm[[group_by_vars[1]]] == groups[group_ind, group_by_vars[1]] &
        preds_wide_perm[[group_by_vars[2]]] == groups[group_ind, group_by_vars[2]] &
        preds_wide_perm[[group_by_vars[3]]] == groups[group_ind, group_by_vars[3]]
      )

      model_scores[preds_inds, model_cols] <- model_scores[preds_inds, base::sample(model_cols)]
    }
  }

  return(model_scores)
}





get_one_perm_each_row_2_models <- function(preds_wide_perm, all_models) {
  model_scores <- preds_wide_perm[, colnames(preds_wide_perm) %in% all_models]

  return(
    apply(model_scores, 1, base::sample) %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(all_models)
  )
}



res <- mean(log_score ~ model + prediction_target + analysis_time_season + region,
  data = preds[preds$before_target_date == "before\ntarget date", ])
temp <- strsplit(names(res), ".", fixed = TRUE) %>%
  unlist() %>%
  matrix(ncol = 4, byrow = TRUE)
res <- cbind(temp, res) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(c("model", "prediction_target", "season", "region", "mean_log_score")) %>%
  `rownames<-`(NULL) %>%
  group_by_(.dots = c("prediction_target", "season", "region")) %>%
  mutate(mean_log_score = as.numeric(mean_log_score),
    rank = rank(-1 * round(as.numeric(mean_log_score), 3)))

res$prediction_target <- factor(res$prediction_target,
  levels = c("Onset Timing", "Peak Timing", "Peak Incidence"))
res$model <- factor(res$model,
  levels = c("KDE", "KCDE", "SARIMA",
    "EW", "CW", "FW-wu",
    "FW-reg-w", "FW-reg-wu", "FW-reg-wui"))

res$region_season <- paste0(res$region, " - ", res$season)

median_per_target_season_region <- res %>%
  group_by(prediction_target, season, region) %>%
  summarize(median_mean_log_score = median(mean_log_score))

res2 <- res %>%
  left_join(median_per_target_season_region, by = c("prediction_target", "season", "region")) %>%
  mutate(diff_from_median = mean_log_score - median_mean_log_score)


res_wide <- res2 %>%
    dplyr::select(-rank, -mean_log_score, -median_mean_log_score) %>%
    tidyr::spread(model, diff_from_median)



perm_test_res <- foreach(i = seq_len(nrow(m_combos))) %dopar% {
  set.seed(60929) # from random.org
  seed_vals <- runif(nrow(m_combos), 1, 10^5)
  set.seed(seed_vals[i])

  model_pair <- m_combos[i, c("m1", "m2")]

  for(j in seq_len(100)) {
    perm_model_scores <-
      lapply(
        1:1000,
        function(k) {
          get_one_perm_each_row_2_models(
            preds_wide_perm = res_wide,
            all_models = model_pair)#,
#            group_by_vars = c("region", "season", "prediction_target"))
        }
      )

    saveRDS(perm_model_scores,
      file = paste0("inst/manuscript/score-permutations/perm_model_scores_season_diff_from_median_",
        model_pair[1],
        "_", model_pair[2],
        "_", j,
        ".rds")
    )
  }
}



res_wide <- res2 %>%
    dplyr::select(-rank, -diff_from_median, -median_mean_log_score) %>%
    tidyr::spread(model, mean_log_score)



perm_test_res <- foreach(i = seq_len(nrow(m_combos))) %dopar% {
  set.seed(60929) # from random.org
  seed_vals <- runif(nrow(m_combos), 1, 10^5)
  set.seed(seed_vals[i])

  model_pair <- m_combos[i, c("m1", "m2")]

  for(j in seq_len(100)) {
    perm_model_scores <-
      lapply(
        1:1000,
        function(k) {
          get_one_perm_each_row_2_models(
            preds_wide_perm = res_wide,
            all_models = model_pair)#,
#            group_by_vars = c("region", "season", "prediction_target"))
        }
      )

    saveRDS(perm_model_scores,
      file = paste0("inst/manuscript/score-permutations/perm_model_scores_season_mean_",
        model_pair[1],
        "_", model_pair[2],
        "_", j,
        ".rds")
    )
  }
}
