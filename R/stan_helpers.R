# R/stan_helpers.R

make_stan_data <- function(df, starting_age) {
  df_aux <- df |>
    dplyr::filter(Age >= starting_age) |>
    dplyr::group_by(Cohort, Sex) |>
    dplyr::mutate(n_Cohort = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(n_Cohort >= 100 - starting_age) |>
    dplyr::mutate(
      Cohort_index = as.numeric(factor(Cohort, levels = sort(unique(Cohort))))
    ) |>
    dplyr::arrange(Cohort, Age)

  if (nrow(df_aux) == 0) return(NULL)

  noise <- runif(nrow(df_aux))
  x <- (df_aux$Age + noise) - starting_age

  list(
    df_aux = df_aux,
    stan_data = list(
      N = length(x),
      x = x,
      D_x = ceiling(df_aux$Dx),
      E_x = df_aux$Ex,
      N_Cohort = max(df_aux$Cohort_index),
      Cohort_idx = df_aux$Cohort_index
    )
  )
}

make_inits <- function(N_Cohort) {
  function() {
    list(
      sigma2 = rep(0.01, N_Cohort),
      b = log(0.1),
      sigma_rw = 0.01,
      trend = 0,
      B_hidden = rep(0, N_Cohort),
      mu_0 = rep(0.001, N_Cohort)
    )
  }
}

fit_model <- function(stan_file, stan_data, control = list(adapt_delta = 0.95, max_treedepth = 12)) {
  rstan::stan(
    file = stan_file,
    data = stan_data,
    init = make_inits(stan_data$N_Cohort),
    control = control,
    save_warmup = FALSE
  )
}

extract_summary <- function(fit, df_aux, sex, country) {
  draws <- rstan::extract(fit)

  data.frame(
    Sex = sex,
    Country = country,
    Cohort = sort(unique(df_aux$Cohort)),
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
}
