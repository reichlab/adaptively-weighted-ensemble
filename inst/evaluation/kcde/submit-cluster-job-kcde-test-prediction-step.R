options(warn = 2)

cores_req <- "1"
mem_req <- "5000"
time_req <- "24:00"
queue_req <- "long"

region_season_combos <- expand.grid(
  region = c("X", paste0("Region-", 1:10)),
  season = paste0(2011:2015, "/", 2012:2016),
  stringsAsFactors = FALSE
)

for(ind in seq_len(nrow(region_season_combos))) {
  region <- region_season_combos$region[ind]
  season <- region_season_combos$season[ind]
  
  output_path <- "/home/er71a/adaptively-weighted-ensemble/inst/evaluation/kcde/cluster-output"
  lsfoutfilename <- "adaptively-weighted-ensemble-kcde-test-predict.out"
  
  case_descriptor <- paste0(region, "-", gsub("/", "-", season))
  filename <- paste0(output_path, "/submit-kcde-test-prediction-step-", case_descriptor, ".sh")
  
  requestCmds <- "#!/bin/bash\n"
  requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
    "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
    "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
    "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
    "#BSUB -W ", time_req, " # run time\n",
    "#BSUB -q ", queue_req, " # which queue we want to run in\n")
  
  cat(requestCmds, file = filename)
  cat("module load R/3.2.2\n", file = filename, append = TRUE)
  cat("module load gsl/1.16\n", file = filename, append = TRUE)
  cat(paste0("R CMD BATCH --vanilla \'--args ",
    region, " ",
    season,
    "\'  /home/er71a/adaptively-weighted-ensemble/inst/evaluation/kcde/kcde-prediction-test.R ",
    output_path, "/output-kcde-test-prediction-step-", case_descriptor, ".Rout"),
    file = filename, append = TRUE)
  
  bsubCmd <- paste0("bsub < ", filename)
  
  system(bsubCmd)
}
