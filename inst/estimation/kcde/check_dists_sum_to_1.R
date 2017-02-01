library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(awes)

loso_preds <- assemble_loso_predictions()

for(prediction_target in c("onset", "peak_week", "peak_inc")) {
  cols_to_examine <- grep(paste0(prediction_target, ".*_log_prob"), colnames(loso_preds), value = TRUE)
  loso_preds[, paste0(prediction_target, "_prob_sum")] <- sapply(
    seq_len(nrow(loso_preds)),
    function(i) {
      sum(exp(loso_preds[i, cols_to_examine]))
    })
}

for(prediction_target in c("ph_1_inc", "ph_2_inc", "ph_3_inc", "ph_4_inc")) {
  cols_to_examine <- grep(paste0(prediction_target, ".*_log_prob"), colnames(loso_preds), value = TRUE)
  loso_preds[, paste0(prediction_target, "_prob_sum")] <- sapply(
    seq_len(nrow(loso_preds)),
    function(i) {
      sum(exp(loso_preds[i, cols_to_examine]))
    })
}

loso_preds$onset_prob_sum
loso_preds$peak_week_prob_sum
loso_preds$peak_inc_prob_sum
loso_preds$ph_1_inc_prob_sum
loso_preds$ph_2_inc_prob_sum
loso_preds$ph_3_inc_prob_sum
loso_preds$ph_4_inc_prob_sum
