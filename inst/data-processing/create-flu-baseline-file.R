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

## fill in baselines for seasons before 2007/2008 with means
mean_region_baseline <- dat %>%
  group_by(region) %>%
  summarize(baseline = round(mean(baseline), 1))
mean_region_baseline <- mean_region_baseline[c(1:2, 4:11, 3), ]

regions <- c("National", paste0("Region", 1:10))
for(first_season_year in 2006:1997) {
  dat <- bind_rows(
    data.frame(
      region = regions,
      season = paste0(first_season_year, "/", first_season_year + 1),
      baseline = mean_region_baseline$baseline
    ),
    dat
  )
}

write.csv(dat, "data-raw/cdc-baselines.csv", quote = FALSE, row.names = FALSE)

## code below here was added and run on 3/1/2017:
flu_onset_baselines <- read.csv("data-raw/cdc-baselines.csv", stringsAsFactors = FALSE)

save(flu_onset_baselines, file = "data/flu_onset_baselines.rdata")
