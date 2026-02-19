# R/02_run_model.R
source(here::here("R", "00_setup.R"))

# Load derived data (or rebuild from HMD if user has it)
data_full <- source(here::here("R", "01_prepare_data.R"))$value

# Settings
sexes <- c("Male", "Female")
country_list <- sort(c(
  "Australia", "Canada", "Denmark", "England & Wales (Total)",
  "Finland", "France (Total)", "Italy", "Japan", "Netherlands",
  "Norway", "Sweden", "United States"
))

starting_ages <- c(80)

# Read Stan model from file
stan_file <- here::here("stan", "gamma_gompertz_rw.stan")

# Helper to prepare one dataset
prepare_stan_data <- function(df, age0) {
  df_aux <- df %>%
    filter(Age >= age0) %>%
    group_by(Cohort, Sex) %>%
    mutate(n_Cohort = n()) %>%
    ungroup() %>%
    filter(n_Cohort >= 100 - age0) %>%
    mutate(Cohort_index = as.numeric(factor(Cohort, levels = sort(unique(Cohort))))) %>%
    arrange(Cohort, Age)

  if (nrow(df_aux) == 0) return(NULL)

  set.seed(1)
  noise <- runif(nrow(df_aux))

  x <- df_aux$Age + noise - age0
  D_x <- ceiling(df_aux$Dx)
  E_x <- df_aux$Ex
  N_Cohort <- max(df_aux$Cohort_index)
  Cohort_idx <- df_aux$Cohort_index

  list(
    stan_data = list(
      N = length(x),
      x = x,
      D_x = D_x,
      E_x = E_x,
      N_Cohort = N_Cohort,
      Cohort_idx = Cohort_idx
    ),
    cohorts = unique(df_aux$Cohort),
    df_aux = df_aux
  )
}

# Reasonable init values (only parameters that exist in Stan)
make_inits <- function(N_Cohort) {
  function() list(
    mu_0 = rep(0.001, N_Cohort),
    sigma2 = rep(0.01, N_Cohort),
    b = log(0.1),
    trend = 0,
    sigma_rw = 0.01,
    B_hidden = rep(0, N_Cohort)
  )
}

dir.create(here::here("outputs"), showWarnings = FALSE)

all_results <- list()

for (age0 in starting_ages) {
  message("=== Starting age: ", age0, " ===")

  stan_model_data <- list()
  sample_post <- NULL

  for (ctry in country_list) {
    for (sx in sexes) {
      message(ctry, " | ", sx, " | age >= ", age0)

      df_sub <- data_full %>%
        filter(Country == ctry, Sex == sx)

      prep <- prepare_stan_data(df_sub, age0)
      if (is.null(prep)) next

      data_list <- prep$stan_data
      N_Cohort <- data_list$N_Cohort

      fit <- rstan::stan(
        file = stan_file,
        data = data_list,
        init = make_inits(N_Cohort),
        control = list(adapt_delta = 0.95, max_treedepth = 12),
        save_warmup = FALSE,
        refresh = 0
      )

      draws <- rstan::extract(fit)

      # Store summary outputs
      aux_fit <- data.frame(
        Sex = sx,
        Country = ctry,
        Cohort = prep$cohorts,
        mu_0 = apply(draws$mu_0, 2, modeest::mlv),
        beta = apply(draws$beta, 2, modeest::mlv),
        B = modeest::mlv(draws$B),
        B_L = HDInterval::hdi(draws$B)[1],
        B_U = HDInterval::hdi(draws$B)[2],
        trend = modeest::mlv(draws$trend),
        trend_L = HDInterval::hdi(draws$trend)[1],
        trend_U = HDInterval::hdi(draws$trend)[2],
        sigma_rw = modeest::mlv(draws$sigma_rw),
        sigma_rw_L = HDInterval::hdi(draws$sigma_rw)[1],
        sigma_rw_U = HDInterval::hdi(draws$sigma_rw)[2],
        sigma2 = apply(draws$sigma2, 2, modeest::mlv)
      )

      sample_post <- rbind(sample_post, aux_fit)

      name_list <- paste(ctry, sx, sep = "_")
      stan_model_data[[name_list]] <- list(
        Sex = sx,
        Country = ctry,
        Cohort = prep$cohorts,
        data = data_list,
        model = fit
      )
    }
  }

  out_file <- here::here("outputs", paste0("model_with_data_", age0, "_outputs.RData"))
  save(data_full, stan_model_data, sample_post, file = out_file)
  message("Saved: ", out_file)

  all_results[[paste0("age", age0)]] <- list(sample_post = sample_post, file = out_file)
}
