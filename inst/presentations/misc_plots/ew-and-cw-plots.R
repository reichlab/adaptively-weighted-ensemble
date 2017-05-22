library(dplyr)
library(ggplot2)

tmp <- data_frame(
    model = rep(c("KDE", "KCDE", "SARIMA"), each=31),
    week = rep(10:40, times=3),
    ew = 1/3,
    cw = rep(c(.2, .5, .3), each=31))

color_palette <- c("#E69F00", "#56B4E9", "#009E73")

ggplot(tmp, aes(x=week, color=model, linetype=model), size=2) +
    geom_line(aes(y=ew), size=1) + theme_bw() + xlab(NULL) +
    ylab("Model Weight") +
    scale_color_manual(breaks = c("KCDE", "KDE", "SARIMA"), values = color_palette) +
    ylim(0,1) 

ggplot(tmp, aes(x=week, color=model, linetype=model)) +
    geom_line(aes(y=cw), size=1) + theme_bw() + xlab(NULL) +
    scale_color_manual(breaks = c("KCDE", "KDE", "SARIMA"), values = color_palette) +
    ylab("Model Weight") +
    ylim(0,1)
