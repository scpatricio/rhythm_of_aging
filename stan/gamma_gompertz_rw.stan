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

  sigma_rw ~ normal(0, 1) T[0, 1]

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
