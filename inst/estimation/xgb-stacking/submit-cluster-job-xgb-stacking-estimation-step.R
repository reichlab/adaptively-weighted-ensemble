options(warn = 2)

all_regions <- c("National", paste0("Region", 1:10))
all_prediction_targets <- c("onset", "peak_week", "peak_inc")
all_explanatory_variables_combos <- c(
  "analysis_time_season_week",
  "analysis_time_season_week-kcde_model_confidence-sarima_model_confidence",
  "analysis_time_season_week-kcde_model_confidence-sarima_model_confidence-weighted_ili"
)
#all_regions <- "National"
#all_prediction_targets <- "onset"

cores_req <- "4"
mem_req <- "15000"
time_req <- "240:00"
#queue_req <- "condo_uma_nicholas_reich"
queue_req <- "condo_grid"

for(region in all_regions) {
  for(prediction_target in all_prediction_targets) {
    for(explanatory_variables in all_explanatory_variables_combos) {
      output_path <- "/home/er71a/adaptively-weighted-ensemble/inst/estimation/xgb-stacking/cluster-output"
      lsfoutfilename <- "adaptively-weighted-ensemble-xgb-stack-est.out"
      
      system_cmd <- paste0("R --vanilla --args ",
        " < inst/estimation/xgb-stacking/xgb-stacking-estimation.R")
      
      case_descriptor <- paste0(
        region,
        "-prediction_target_", prediction_target,
        "_", explanatory_variables
      )
      filename <- paste0(output_path, "/submit-xgb-stacking-estimation-step-", case_descriptor, ".sh")
      
      requestCmds <- "#!/bin/bash\n"
      requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
        "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
        "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
        "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
        "#BSUB -W ", time_req, " # run time\n",
        "#BSUB -q ", queue_req, " # which queue we want to run in\n")
      
      cat(requestCmds, file = filename)
      cat("module load openmpi/2.0.1\n", file = filename, append = TRUE)
      cat("module load R/3.2.2\n", file = filename, append = TRUE)
      cat(paste0("R CMD BATCH --vanilla \'--args ",
        region, " ",
        prediction_target, " ",
        explanatory_variables, " ",
        cores_req,
        "\'  /home/er71a/adaptively-weighted-ensemble/inst/estimation/xgb-stacking/xgb-stacking-estimation.R ",
        output_path, "/output-xgb-stacking-estimation-", case_descriptor, ".Rout"),
        file = filename, append = TRUE)
      
      bsubCmd <- paste0("bsub < ", filename)
      
      system(bsubCmd)
    } # explanatory_variables
  } # prediction_target
} # region
