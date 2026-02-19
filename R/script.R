###############################################################################
# Project: Gamma–Gompertz + Random Walk (with drift) decomposition
# Purpose: Fit cohort-specific gamma–Gompertz hazards with a latent RW-on-log(b)
# Output: Saved fitted stan objects + posterior summaries by age threshold
###############################################################################

# ----------------------------
# 0) Session setup
# ----------------------------

setwd("~/Desktop/bayesian/Mixture_model/Period/GG_spline/life_exp/")
set.seed(1)

rm(list = ls(all = TRUE))

# Core libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)

library(ggplot2)
library(splines)

# Bayesian modeling
library(rstan)

# Parallel backend
library(doParallel)

# ----------------------------
# 1) Parallel setup (Stan / R)
# ----------------------------

totalCores = detectCores()

# Leave one core free to avoid overloading the machine
cluster = makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ----------------------------
# 2) Stan model specification
# ----------------------------
# Model:
#   - Cohort-specific baseline: mu_0[t]
#   - Cohort-specific heterogeneity: sigma2[t]
#   - Cohort-specific Gompertz slope: beta[t] = exp(b + B_hidden[t] - mean(B_hidden))
#   - Latent cohort process: B_hidden follows a random walk with drift (Laplace innovations)
#   - Deaths modeled as Poisson(D_x | lambda(x,t) * exposure)
#
# NOTE: mu_0 and sigma2 are constrained to [0,1] here. This is your original choice.
# NOTE: sigma_rw is constrained to [0,1] via T[0,1] in the prior. Original choice retained.

model_string = "
  data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0> D_x[N];
  vector[N] E_x;

  int<lower=0> N_Cohort;
  int<lower=1, upper=N_Cohort> Cohort_idx[N];
}
parameters {
  vector<lower=0, upper=1>[N_Cohort] mu_0;
  vector<lower=0, upper=1>[N_Cohort] sigma2;

  real b;
  vector[N_Cohort] B_hidden;
  real trend;
  real<lower=0> sigma_rw;
}
transformed parameters {
  vector[N] lambda;
  vector[N_Cohort] beta = exp(b + B_hidden - mean(B_hidden));
  real B = exp(b);

  for (i in 1:N) {
    int t = Cohort_idx[i];
    lambda[i] = mu_0[t] * exp(beta[t] * x[i]) /
      (1 + sigma2[t] * mu_0[t] * (exp(beta[t] * x[i]) - 1) / beta[t]);
  }
}
model {
  mu_0 ~ normal(0, 1) T[0, 1];
  sigma2 ~ gamma(1, 0.5);
  b ~ normal(0, 2);
  trend ~ normal(0, 2);

  // Random walk prior (Laplace innovations)
  B_hidden[1] ~ double_exponential(trend, sigma_rw);
  B_hidden[2:N_Cohort] ~ double_exponential(B_hidden[1:(N_Cohort-1)] + trend, sigma_rw);

  sigma_rw ~ normal(0, 1) T[0, 1];

  D_x ~ poisson(lambda .* E_x);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(D_x[i] | lambda[i] * E_x[i]);
    y_rep[i] = poisson_rng(lambda[i] * E_x[i]);
  }
}
"

# ----------------------------
# 3) Data (HMD cohort-format)
# ----------------------------

data_full = read.csv2("~/Downloads/rhythm_of_aging-main/data/PNAS_cohort_data_HMD.csv")

# ----------------------------
# 4) Country/sex selection and loop setup
# ----------------------------

sample_post = NULL
stan_model_data = list()

# Countries to include in this run
countries_include = c('Denmark', 'Finland', 'France (Total)', 'Italy',
                      'Netherlands', 'England & Wales (Total)', 'Norway', 
                      'Sweden', 'Australia', "Canada", 'United States',
                      'Japan')

country_list = sort(country_list)

# Sexes
sexes = c("Male", "Female")

# NOTE: The next few assignments appear to be exploratory/debug remnants,
# but they do not affect the final loops because they are overwritten later.
country_loop = sample(country_list, 1)
sex_loop = "Male"

# set initical age
age_0 = 80

# ----------------------------
# 5) Main estimation loop
#    Fit model for multiple starting ages, countries, and sexes
# ----------------------------

stan_model_data = list()
sample_post = NULL

for (country_loop in countries_include) {
  for (sex_loop in sexes) {
    
    print(paste(country_loop, sex_loop, "from age", age_0))
    
    # ----------------------------
    # 5.1) Subset and filter data
    # ----------------------------
    # - Keep ages >= age_0
    # - Require minimum age coverage per cohort: n_Cohort >= (100 - age_0)
    #   (i.e., ~complete coverage up to age 100)
    # - Construct ordered cohort index for Stan
    df_aux = data_full %>%
      filter(
        Country == country_loop,
        Sex == sex_loop,
        Age >= age_0
      ) %>%
      group_by(Cohort, Sex) %>%
      mutate(n_age_groups = n()) %>%
      ungroup() %>%
      filter(n_age_groups >= 100 - age_0) %>%
      mutate(
        Cohort_index = as.numeric(factor(Cohort, levels = sort(unique(Cohort))))
      ) %>%
      arrange(Cohort, Age)
    
    if (nrow(df_aux) == 0) next
    
    # ----------------------------
    # 5.2) Prepare Stan inputs
    # ----------------------------
    
    # Add small uniform noise to age to break ties and improve smoother plotting
    # (and to avoid identical x within cohort if needed)
    set.seed(1)
    noise = runif(nrow(df_aux))
    
    # x is re-centered to start at 0 at age_0
    x = df_aux$Age + noise - age_0
    
    D_x = ceiling(df_aux$Dx)
    E_x = df_aux$Ex
    
    N_Cohort = max(df_aux$Cohort_index)
    Cohort_idx = df_aux$Cohort_index
    
    # ----------------------------
    # 5.3) Quick visualization (optional diagnostic)
    # ----------------------------
    # df_aux %>%
    #   mutate(Age = Age + noise) %>%
    #   filter(Dx > 0) %>%
    #   ggplot(aes(x = Age, y = mx, group = Cohort, color = Cohort)) +
    #   geom_point(alpha = 0.3) +
    #   geom_smooth(se = FALSE) +
    #   scale_color_viridis_c() +
    #   scale_y_log10() +
    #   theme_minimal() +
    #   guides(color = "none")
    
    # Stan data list
    data_list = list(
      N = length(x),
      x = x,
      D_x = D_x,
      E_x = E_x,
      N_Cohort = N_Cohort,
      Cohort_idx = Cohort_idx
      )
    
    # ----------------------------
    # 5.4) Initial values
    # ----------------------------
    init_list = function() {
      list(
        sigma2 = rep(0.01, N_Cohort),
        b = log(0.1),
        sigma_rw = 0.01,
        trend = 0,
        B_hidden = rep(0, N_Cohort),
        mu_0 = rep(0.001, N_Cohort)
      )
    }
    
    # ----------------------------
    # 5.5) Fit model in Stan
    # ----------------------------
    warmup = 4000
    
    fit = stan(
      model_code = model_string,
      data = data_list,
      chains = 4,
      iter = warmup + 2000,
      warmup = warmup,
      init = init_list,
      # refresh = 0,
      save_warmup = FALSE
    )
    
    # ----------------------------
    # 5.6) Posterior predictive checks
    # ----------------------------
    
    draws = rstan::extract(fit)
    y_rep = draws$y_rep  # [iterations, N]
    
    bayesplot::ppc_dens_overlay(
      y = D_x,
      yrep = y_rep[sample(1:nrow(y_rep), 100), ],
      main = "Posterior Predictive Check: Observed vs. Simulated Death Counts"
    ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # ----------------------------
    # 5.7) Key parameter summaries and diagnostics
    # ----------------------------
    
    parameters = c("B", "trend", "sigma_rw")
    to_print = summary(fit, pars = parameters)$summary
    print(to_print)
    
    plot(fit, pars = "beta")
    plot(fit, pars = "mu_0")
    plot(fit, pars = "sigma2")
    
    # Latent RW diagnostics
    B_hidden = apply(draws$B_hidden, 2, median)
    plot(B_hidden, type = "l")
    plot(diff(B_hidden), type = "l")
    
    beta = apply(draws$beta, 2, median)
    plot(beta, type = "l")
    
    mu_0 = apply(draws$mu_0, 2, median)
    plot((mu_0), type = "l")
    
    to_print = summary(fit, pars = parameters)$summary
    print(to_print)
    
    traceplot(fit, pars = parameters, inc_warmup = FALSE)
    stan_dens(fit, pars = parameters, inc_warmup = FALSE)
    
    # HDI of a derived quantity (kept as provided)
    HDInterval::hdi(apply(draws$beta, 1, function(x) exp(modeest::mlv(log(x)))))
    
    # ----------------------------
    # 5.8) Extract posterior summaries for storage
    # ----------------------------
    
    fit_params = rstan::extract(fit)
    
    aux_fit = data.frame(
      Sex = sex_loop,
      Country = country_loop,
      Cohort = unique(df_aux$Cohort),
      
      mu_0 = apply(fit_params$mu_0, 2, modeest::mlv),
      beta = apply(fit_params$beta, 2, modeest::mlv),
      
      B = modeest::mlv(fit_params$B),
      B_L = HDInterval::hdi(fit_params$B)[1],
      B_U = HDInterval::hdi(fit_params$B)[2],
      
      trend = modeest::mlv(fit_params$trend),
      trend_L = HDInterval::hdi(fit_params$trend)[1],
      trend_U = HDInterval::hdi(fit_params$trend)[2],
      
      sigma_rw = modeest::mlv(fit_params$sigma_rw),
      sigma_rw_L = HDInterval::hdi(fit_params$sigma_rw)[1],
      sigma_rw_U = HDInterval::hdi(fit_params$sigma_rw)[2],
      
      sigma2 = apply(fit_params$sigma2, 2, modeest::mlv)
    )
    
    sample_post = rbind(sample_post, aux_fit)
    
    # ----------------------------
    # 5.9) Save fit object and data used
    # ----------------------------
    
    name_list = paste(country_loop, sex_loop, sep = "_")
    
    list_fit = list(
      Sex = sex_loop,
      Country = country_loop,
      Cohort = unique(df_aux$Cohort),
      data = data_list,
      model = fit
    )
    
    stan_model_data[[name_list]] = list_fit
  }
}

# ----------------------------
# 6) Save outputs per starting age
# ----------------------------

file_name = paste0("model_with_data_", age_0, "_new_model_outputs.RData")
save(data_full, stan_model_data, sample_post, file = file_name)
