## download and clean regional US flu data
## October 5 2016 - Nicholas Reich - created
## October 7 2016 - Evan Ray - calculate time using MMWRweek package, fix bug in
##   computation of season week
## October 7 2016 - Nicholas Reich - merged US and Regional data into one file.

library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(MMWRweek)
library(cdcfluview)

regionflu <- get_flu_data("hhs", sub_region=1:10, data_source="ilinet", years=1997:2016)
usflu <- get_flu_data("national", sub_region=NA, data_source="ilinet", years=1997:2016)

## make AGE cols in usflu integer data type
cols <- matches('^AGE', vars=colnames(usflu))
usflu[,cols] <- sapply(usflu[,cols], as.integer)
cols <- matches('^AGE', vars=colnames(regionflu))
regionflu[,cols] <- sapply(regionflu[,cols], as.integer)

data <- bind_rows(regionflu, usflu)

data <- transmute(data,
  region.type = `REGION TYPE`,
  region = REGION,
  year = YEAR,
  week = WEEK,
  time = as.POSIXct(MMWRweek2Date(YEAR, WEEK)),
  weighted_ili = as.numeric(`% WEIGHTED ILI`))

## set zeroes to NAs
data[which(data$weighted_ili==0),"weighted_ili"] <- NA

## Add time_index column: the number of days since some origin date (1970-1-1 in this case).
## The origin is arbitrary.
data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))

## Season column: for example, weeks of 2010 up through and including week 30 get season 2009/2010;
## weeks after week 30 get season 2010/2011
## I am not sure where I got that the season as defined as starting on MMWR week 30 from...
data$season <- ifelse(
  data$week <= 30,
  paste0(data$year - 1, "/", data$year),
  paste0(data$year, "/", data$year + 1)
)

## Season week column: week number within season
## weeks after week 30 get season_week = week - 30
## weeks before week 30 get season_week = week + (number of weeks in previous year) - 30
## This computation relies on the start_date function in package MMWRweek,
## which is not exported from that package's namespace!!!
data$season_week <- ifelse(
  data$week <= 30,
  data$week + MMWRweek(MMWRweek:::start_date(data$year) - 1)$MMWRweek - 30,
  data$week - 30
)

write.csv(data, file = "data-raw/allflu-cleaned.csv")
