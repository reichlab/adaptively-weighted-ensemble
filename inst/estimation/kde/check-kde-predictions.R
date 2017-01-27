library(tidyr)
library(ggplot2)
library(dplyr)

region_strings <- c("National", paste0("Region", 1:10))

pdf("inst/estimation/kde/check-kde-predictions.pdf", width=10)
for(reg in region_strings) {
    analysis_time_seasons <- paste0(1997:2010, "-", 1998:2011)
    ## read in first file to set file standard
    fname <- paste0("inst/estimation/loso-predictions/kde-", reg, "-", 
                    analysis_time_seasons[1], "-loso-predictions.rds")
    tmp <- readRDS(fname)
    
    for(i in 2:length(analysis_time_seasons)) {
        fname <- paste0("inst/estimation/loso-predictions/kde-", reg, "-", 
                        analysis_time_seasons[i], "-loso-predictions.rds")
        tmp_to_add <- readRDS(fname)
        tmp <- bind_rows(tmp, tmp_to_add)
    }
    
    tmp_long <- as_data_frame(tmp) %>% 
        select(model, analysis_time_season, analysis_time_season_week,
               starts_with("prediction_week_ph"), 
               starts_with("analysis_time"), contains("log_score"),
               -contains("competition")) %>%
        gather(key=metric, value=log_score, -c(model, 
                                               starts_with("prediction_week_ph"), 
                                               starts_with("analysis_time"))) %>%
        ## exclude pandemic season
        filter(analysis_time_season != "2009/2010")
    
    p <- ggplot(tmp_long, aes(x=analysis_time_season_week, y=log_score)) +
        geom_line(aes(color=analysis_time_season)) +
        facet_grid(.~metric) +
        geom_smooth(se=FALSE, color="black") +
        ylim(-10, 0)
    print(p)
}
dev.off()