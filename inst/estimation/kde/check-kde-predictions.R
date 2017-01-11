library(tidyr)
library(ggplot2)
library(dplyr)

region_strings <- c("National", paste0("Region", 1:10))

pdf("inst/estimation/kde/check-kde-predictions.pdf", width=10)
for(reg in region_strings) {
    fname <- paste0("inst/estimation/loso-predictions/kde-", reg, "-loso-predictions.rds")
    tmp <- readRDS(fname)
    tmp <- as_data_frame(tmp) %>% 
        gather(key=metric, value=log_score, -c(model, 
                                               starts_with("prediction_week_ph"), 
                                               starts_with("analysis_time"))) %>%
        ## exclude pandemic season
        filter(analysis_time_season != "2009/2010")
    
    p <- ggplot(tmp, aes(x=analysis_time_season_week, y=log_score)) +
        geom_line(aes(color=analysis_time_season)) +
        facet_grid(.~factor(metric)) +
        geom_smooth(se=FALSE, color="black") +
        ylim(-10, 0)
    print(p)
}
dev.off()