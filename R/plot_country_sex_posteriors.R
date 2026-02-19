# ============================================================
# Plot posterior summaries (b, drift, sigma_rw) across countries
# ============================================================
# This script:
#   1) Loads saved Stan fits (stan_model_data) from an .RData file
#   2) Extracts posterior summaries for each country-sex series
#   3) Builds a 3-panel figure: b, drift, sigma_rw with 95% HDIs
#
# Notes:
#   - Assumes 'stan_model_data' exists in the loaded .RData file.
#   - Uses posterior mode (via modeest::mlv) and HDInterval::hdi.
#   - Intended for reproducible plotting for the PNAS manuscript.
# ============================================================

set.seed(1)
rm(list = ls(all = TRUE))

# -----------------------------
# User inputs (edit as needed)
# -----------------------------
model_rdata_path <- "model_outputs_with_data.RData"  # path to saved fits

# -----------------------------
# Packages
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(tidyr)

  # Needed for extraction and summaries
  library(rstan)
  library(modeest)
  library(HDInterval)
})

# -----------------------------
# Load fitted models
# -----------------------------
load(model_rdata_path)
# Expecting: stan_model_data (list of lists with $model, $Sex, $Country)
full_stan_model_data <- stan_model_data

# -----------------------------
# Helper functions
# -----------------------------

# Extract posterior summaries from one fitted model object
extract_series_summaries <- function(fit_entry) {
  posterior <- rstan::extract(fit_entry$model)

  # This reproduces your original logic:
  # For each iteration: compute a "typical" beta value by taking
  # exp(mode(log(beta))) and then take the mode across iterations.
  beta_iter_summary <- apply(
    posterior$beta,
    1,
    function(x) exp(modeest::mlv(log(x)))
  )

  data.frame(
    Sex      = fit_entry$Sex,
    Country  = fit_entry$Country,

    B        = modeest::mlv(beta_iter_summary),
    B_L      = HDInterval::hdi(beta_iter_summary)[1],
    B_U      = HDInterval::hdi(beta_iter_summary)[2],

    trend    = modeest::mlv(posterior$trend),
    trend_L  = HDInterval::hdi(posterior$trend)[1],
    trend_U  = HDInterval::hdi(posterior$trend)[2],

    sigma_rw   = modeest::mlv(posterior$sigma_rw),
    sigma_rw_L = HDInterval::hdi(posterior$sigma_rw)[1],
    sigma_rw_U = HDInterval::hdi(posterior$sigma_rw)[2]
  )
}

# Extract ggplot legend
extract_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

# -----------------------------
# Build summary dataframe
# -----------------------------
B_df_list <- lapply(full_stan_model_data, extract_series_summaries)
B_df <- do.call(rbind, B_df_list)
rownames(B_df) <- NULL

# Clean country labels
B_df$Country <- sub(" \\(Total\\)", "", B_df$Country)

# Order countries alphabetically (reversed for top-to-bottom readability)
B_df$Country <- factor(B_df$Country, levels = rev(sort(unique(B_df$Country))))

# -----------------------------
# Plot styling parameters
# -----------------------------
subtitle_size <- 12
country_size  <- 10
x_axis_size   <- 10

sex_colors <- c("Male" = "steelblue", "Female" = "firebrick")

# -----------------------------
# Panel 1: Rate of aging (b)
# -----------------------------
p1 <- ggplot(B_df, aes(x = B, y = Country, color = Sex)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = B_L, xmax = B_U),
    width = 0.2,
    position = position_dodge(width = 0.5),
    size = 0.7
  ) +
  scale_color_manual(values = sex_colors) +
  labs(x = "", y = "", subtitle = expression(italic("Rate of aging ("*b*")"))) +
  theme_minimal() +
  theme(
    axis.text.x    = element_text(hjust = 0.5, size = x_axis_size),
    axis.text.y    = element_text(size = country_size),
    legend.position = "bottom",
    plot.subtitle  = element_text(size = subtitle_size, face = "italic")
  )

# -----------------------------
# Panel 2: Drift (beta)
# -----------------------------
p2 <- ggplot(B_df, aes(x = trend, y = Country, color = Sex)) +
  geom_vline(xintercept = 0, col = "black", size = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = trend_L, xmax = trend_U),
    width = 0.2,
    position = position_dodge(width = 0.5),
    size = 0.7
  ) +
  scale_color_manual(values = sex_colors) +
  labs(x = "", y = "", subtitle = expression(italic("Drift ("*beta*")"))) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(hjust = 0.5, size = x_axis_size),
    plot.subtitle = element_text(size = subtitle_size, face = "italic")
  )

# -----------------------------
# Panel 3: Period volatility (sigma_rw)
# -----------------------------
p3 <- ggplot(B_df, aes(x = sigma_rw, y = Country, color = Sex)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(xmin = sigma_rw_L, xmax = sigma_rw_U),
    width = 0.2,
    position = position_dodge(width = 0.5),
    size = 0.7
  ) +
  scale_x_continuous(breaks = c(0.025, 0.075, 0.125)) +
  scale_color_manual(values = sex_colors) +
  labs(x = "", y = "", subtitle = expression(italic("Period effect (" * sigma[rw] * ")"))) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(hjust = 0.5, size = x_axis_size),
    plot.subtitle = element_text(size = subtitle_size, face = "italic")
  )

# -----------------------------
# Combine panels with one legend
# -----------------------------
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none", axis.text.y = element_blank())
p3_noleg <- p3 + theme(legend.position = "none", axis.text.y = element_blank())

combined_panels <- (p1_noleg | p2_noleg | p3_noleg) +
  plot_layout(widths = c(1, 1, 1))

shared_legend <- extract_legend(
  p1 + theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10),
    legend.margin    = margin(t = -40),
    legend.box.margin = margin(t = 10)
  )
)

final_plot <- cowplot::plot_grid(
  combined_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
) +
  theme(plot.margin = margin(t = -8, r = -4, b = -10, l = -23))

# -----------------------------
# Display
# -----------------------------
final_plot
