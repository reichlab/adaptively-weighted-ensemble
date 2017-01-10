## make a single prediction from KDE/GAM model
## Nicholas Reich
## November 2016

library(KoTflu20162017)

source("inst/estimation/kde/kde-utils.R")
source("inst/estimation/kde/kde-prediction-utils.R")
source("inst/estimation/kde/make_one_kde_prediction_file.R")

n_sims <- 10000

## for 2016-2017 season, first predictions due on 11/7/2016 (EW45 == SW15)
## using data posted on 11/4/2016 that includes up to EW43 == SW13
## last predictions due on 5/15/2017 (EW20 == SW 42)
## last predictions use data through EW40 == SW18
## first_analysis_time_season_week could be set to 13, but padding at front end

for(season_week in 13:42) {
    make_one_kde_prediction_file(path = "inst/estimation/kde/current-fits/",
                                 season = "2016/2017",
                                 season_week = season_week,
                                 n_sim = n_sims)
}
