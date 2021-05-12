// congeneric model (heteroskedastic) with log-normal loadings
data {
  int<lower = 2> Np;
  int<lower = 1> Ni;
  int<lower = 2> N;
  real<lower = 0> alpha_scale;
  real<lower = 0> beta_scale;
  real lambda_median;
  real<lower = 0> lambda_scale;
  int<lower = 1, upper = Np> resp_id[N];
  int<lower = 1, upper = Ni> item_id[N];
  vector[N] y;
  int ret_yhat;
  int ret_ll;
}
transformed data {
  int Ny = 0;
  int Nll = 0;

  if (ret_yhat == 1) Ny = N;
  if (ret_ll == 1) Nll = N;
}
parameters {
  real alpha;
  real<lower = 0> sigma_beta;
  vector[Ni] beta;
  vector<lower = 0>[Ni] lambda;
  real ln_alpha;
  real<lower = 0> ln_sigma_beta;
  vector[Ni] ln_sigma2;
  vector<lower = 0>[Ni] ln_lambda;
  vector[2] Theta[Np];
  real<lower = 0, upper = 1> r_tr;
}
model {
  matrix[2, 2] R;

  R[1, 1] = 1;
  R[2, 2] = 1;
  R[1, 2] = r_tr * 2 - 1;
  R[2, 1] = r_tr * 2 - 1;

  r_tr ~ beta(2, 2);

  alpha ~ normal(0, alpha_scale);
  sigma_beta ~ normal(0, beta_scale);
  beta ~ normal(alpha, sigma_beta);

  ln_alpha ~ normal(0, 5);
  ln_sigma_beta ~ normal(0, 5);
  ln_sigma2 ~ normal(ln_alpha, ln_sigma_beta);

  lambda ~ lognormal(lambda_median, lambda_scale);
  ln_lambda ~ lognormal(lambda_median, lambda_scale);

  Theta ~ multi_normal(rep_vector(0, 2), R);

  {
    vector[N] sigma_y;
    vector[N] mu;
    for (i in 1:N) {
      sigma_y[i] = sqrt(exp(ln_sigma2[item_id[i]] + ln_lambda[item_id[i]] * Theta[resp_id[i], 2]));
      mu[i] = beta[item_id[i]] + lambda[item_id[i]] * Theta[resp_id[i], 1];
    }
    y ~ normal(mu, sigma_y);
  }
}
generated quantities {
  vector[Nll] log_lik;
  vector[Ny] yhat;
  vector[Ni] sigma = sqrt(exp(ln_sigma2));
  real r = r_tr * 2 - 1;

  {
    real sigma_y;
    real mu;
    for (i in 1:max(Nll, Ny)) {
      sigma_y = sqrt(exp(ln_sigma2[item_id[i]] + ln_lambda[item_id[i]] * Theta[resp_id[i], 2]));
      mu = beta[item_id[i]] + lambda[item_id[i]] * Theta[resp_id[i], 1];
      if (Nll > 0) log_lik[i] = normal_lpdf(y[i] | mu, sigma_y);
      if (Ny > 0) yhat[i] = normal_rng(mu, sigma_y);
    }
  }
}
