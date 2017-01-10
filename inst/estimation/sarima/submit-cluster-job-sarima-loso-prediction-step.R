options(warn = 2, error = recover)

cores_req <- "1"
mem_req <- "5000"
time_req <- "12:00"
queue_req <- "long"

region_season_combos <- expand.grid(
  region = c("X", paste0("Region ", 1:10)),
  season = paste0(1997:2010, "/", 1998:2011),
  stringsAsFactors = FALSE
)

for(ind in seq_len(nrow(region_season_combos))) {
  region <- region_season_combos$region[ind]
  season <- region_season_combos$season[ind]
  
  save_path <- "/home/er71a/adaptively-weighted-ensemble/inst/estimation/sarima/fits"
  output_path <- "/home/er71a/adaptively-weighted-ensemble/inst/estimation/sarima/cluster-estimation-output"
  lsfoutfilename <- "adaptively-weighted-ensemble-sarima-predict.out"
  
  requestCmds <- "#!/bin/bash\n"
  requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
    "#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
    "#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
    "#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
    "#BSUB -W ", time_req, " # run time\n",
    "#BSUB -q ", queue_req, " # which queue we want to run in\n")
  
  cat(requestCmds, file = filename)
  cat("module load R/3.2.2\n", file = filename, append = TRUE)
  cat(paste0("R CMD BATCH --vanilla \'--args ",
    region, " ",
    season,
    "\'  /home/er71a/adaptively-weighted-ensemble/inst/estimation/sarima/sarima-prediction-loso.R ",
    output_path, "/output-sarima-estimation-step-", case_descriptor, ".Rout"),
    file = filename, append = TRUE)
  
  bsubCmd <- paste0("bsub < ", filename)
  
  system(bsubCmd)
}
