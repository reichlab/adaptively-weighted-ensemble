## create long flu baseline dataset
## Nicholas Reich
## October 2016

library(tidyr)
library(dplyr)

dat <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/wILI_Baseline.csv")
dat <- as_data_frame(dat) %>%
    gather(key=year, value=baseline, -X) %>%
    transmute(region = as.character(X),
              season = substr(year, start=2, stop=10), ## removes X at beginning of season name
              baseline=baseline)
dat$season <- gsub("\\.", "/", dat$season) ## replaces . with / in season string

write.csv(dat, "data-raw/cdc-baselines.csv", quote = FALSE, row.names = FALSE)