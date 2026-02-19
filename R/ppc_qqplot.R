# ============================================================
# Posterior predictive QQ diagnostics (death counts)
# ============================================================
# This script:
#   1) Loads saved Stan fits (stan_model_data) from an .RData file
#   2) Builds a quantile–quantile diagnostic comparing:
#        - empirical quantiles from observed deaths (D_x)
#        - posterior predictive quantiles from y_rep
#      for each country-sex series
#   3) Plots QQ curves faceted by sex
#
# Notes:
#   - Assumes 'stan_model_data' exists in the .RData file and that each
#     list element has: $model (Stan fit), $data$D_x (observed deaths),
#     plus metadata $Sex and $Country.
#   - The QQ calculation reproduces the original logic exactly.
# ============================================================

rm(list = ls(all = TRUE))

# -----------------------------
# User inputs (edit as needed)
# -----------------------------
model_rdata_path <- "model_outputs_with_data.RData"

# -----------------------------
# Packages
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Needed to extract posterior draws
  library(rstan)
})

# -----------------------------
# Load fitted models
# -----------------------------
load(model_rdata_path)
full_stan_model_data <- stan_model_data

# -----------------------------
# Helper function: QQ data for one country-sex series
# -----------------------------
compute_ppc_qq <- function(fit_entry) {
  posterior <- rstan::extract(fit_entry$model)

  # Posterior predictive replicated deaths:
  # y_rep dimensions are typically [iterations, N]
  y_rep <- posterior$y_rep

  # Observed deaths used in the likelihood
  y_obs <- fit_entry$data$D_x

  # Quantile grid: same as original
  qt <- seq(0.01, 0.99, by = 0.01)
  qt <- qt[-1]

  # Empirical quantiles of observed deaths (D_x)
  dx_obs <- quantile(y_obs, probs = qt)

  # ------------------------------------------------------------
  # IMPORTANT:
  # The following reproduces your original logic exactly:
  #
  #   q_post <- mean(y_rep <= dx)
  #
  # Since y_rep is a matrix, this returns the proportion of all
  # replicated values (across all iterations and observations)
  # that are <= dx. This is *not* the same as computing a per-i
  # posterior predictive quantile. Keeping as-is by request.
  # ------------------------------------------------------------
  qq_mat <- sapply(dx_obs, function(dx) {
    q_emp  <- mean(y_obs <= dx)
    q_post <- mean(y_rep <= dx)
    c(dx, q_emp, q_post)
  })

  qq_df <- as.data.frame(t(qq_mat))
  names(qq_df) <- c("qt", "observed", "posterior")

  qq_df$Sex     <- fit_entry$Sex
  qq_df$Country <- fit_entry$Country

  qq_df
}

# -----------------------------
# Build combined QQ dataset
# -----------------------------
qq_list <- lapply(full_stan_model_data, compute_ppc_qq)
qq_df <- do.call(rbind, qq_list)
rownames(qq_df) <- NULL

# Clean country labels (remove "(Total)")
qq_df$Country <- sub(" \\(Total\\)", "", qq_df$Country)

# Order countries consistently (reversed alphabetical)
qq_df$Country <- factor(qq_df$Country, levels = rev(sort(unique(qq_df$Country))))

# -----------------------------
# Plot: Posterior vs observed quantiles
# -----------------------------
subtitle_size <- 12
legend_size   <- 11
axis_size     <- 10

qq_plot <- ggplot(qq_df, aes(x = observed, y = posterior, color = Country)) +
  geom_abline(slope = 1, linetype = "dashed", size = 1, alpha = 0.2) +
  geom_line(alpha = 0.5, size = 0.8) +
  facet_wrap(~Sex, nrow = 2) +
  labs(
    x = "Observed Quantile",
    y = "Posterior Quantile",
    color = "",
    subtitle = "QQ-Plot"
  ) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(hjust = 0.5, size = axis_size),
    axis.text.y   = element_text(size = axis_size),
    legend.text   = element_text(size = legend_size),
    legend.position = "right",
    plot.subtitle = element_text(size = subtitle_size, face = "italic")
  )

print(qq_plot)
