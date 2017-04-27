library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(awes)

obs_peaks_training <- sapply(paste0(1997:2010, "/", 1998:2011),
  function(season) {
    get_observed_seasonal_quantities(
      data = flu_data %>%
        filter(region == "X"),
      season = season,
      first_CDC_season_week = 10,
      last_CDC_season_week = 42,
      onset_baseline = 2.1,
      incidence_var = "weighted_ili",
      incidence_bins = data.frame(
        lower = c(0, seq(from = 0.05, to = 12.95, by = 0.1)),
        upper = c(seq(from = 0.05, to = 12.95, by = 0.1), Inf)),
      incidence_bin_names = as.character(seq(from = 0, to = 13, by = 0.1))
    )$observed_peak_inc
  }
) %>%
  unlist() %>%
  unname()

peaks_dens_fit <- density(obs_peaks_training)
bw <- peaks_dens_fit$bw

x_grid_size <- 101
x_grid_range <- c(min(obs_peaks_training) - 4 * bw,
  max(obs_peaks_training) + 4 * bw)
x_grid <- seq(from = x_grid_range[1],
  to = x_grid_range[2],
  length = x_grid_size)

# contributions from individual seasons
plot_df <- rbind.fill(
  lapply(obs_peaks_training,
    function(peak_val) {
      data.frame(
        x = x_grid,
        y = dnorm(x_grid, mean = peak_val, sd = bw) / length(obs_peaks_training),
        group = peak_val,
        type = "Contribution"
      )
    })
  )

# combined density estimate
combined_vals <- plot_df %>%
  group_by(x) %>%
  summarize(y = sum(y)) %>%
  mutate(group = "Combined",
    type = "Combined")

# put them together
plot_df <- rbind(
  plot_df,
  combined_vals
  ) %>%
  mutate(group = factor(as.character(group)))

peaks_df <- data.frame(
  x = obs_peaks_training,
  y = 0
)

p <- ggplot(plot_df) +
  geom_line(aes(x = x, y = y, group = group, colour = type, linetype = type)) +
  geom_point(aes(x = x, y = y), data = peaks_df) +
  scale_colour_manual("Density",
    breaks = c("Contribution", "Combined"),
    values = c("grey", "black")) +
  scale_linetype_manual("Density",
    breaks = c("Contribution", "Combined"),
    values = c(2, 1)) +
  xlab("Peak Incidence") +
  ylab("Density") +
  theme_bw()
print(p)
