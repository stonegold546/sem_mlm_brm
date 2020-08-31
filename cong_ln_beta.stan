// congeneric with log-normal loadings
// Distribution is beta
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
  real scaler;
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
  vector[Np] theta_p;
  vector<lower = 0>[Ni] inv_prec;
  real prec_mean;
  real<lower = 0> prec_scale;
}
transformed parameters {
  vector[Ni] prec = 1 ./ inv_prec;
}
model {
  alpha ~ normal(0, alpha_scale);
  sigma_beta ~ normal(0, beta_scale);
  beta ~ normal(alpha, sigma_beta);

  lambda ~ lognormal(lambda_median, lambda_scale);
  theta_p ~ std_normal();

  inv_prec ~ lognormal(prec_mean, prec_scale);
  prec_mean ~ std_normal();
  prec_scale ~ std_normal();

  {
    vector[N] prob = inv_logit(beta[item_id] + lambda[item_id] .* theta_p[resp_id]);
    vector[N] prec_y = prec[item_id];
    y ~ beta(prob .* prec_y, (1 - prob) .* prec_y);
  }
}
generated quantities {
  matrix[Ni, 2] beta_pars;
  vector[Ni] i_means;
  vector[Ni] sigma;
  vector[Nll] log_lik;
  vector[Ny] yhat;

  {
    real p;
    for (i in 1:max(Nll, Ny)) {
      p = inv_logit(beta[item_id[i]] + lambda[item_id[i]] * theta_p[resp_id[i]]);
      if (Nll > 0) log_lik[i] = beta_lpdf(y[i] | p * prec[item_id[i]], (1 - p) * prec[item_id[i]]);
      if (Ny > 0) yhat[i] = beta_rng(p * prec[item_id[i]], (1 - p) * prec[item_id[i]]) * scaler;
    }
  }

  {
    vector[Ni] p = inv_logit(beta);
    i_means = p * scaler;
    beta_pars[, 1] = p .* prec;
    beta_pars[, 2] = (1 - p) .* prec;
    sigma = sqrt(p .* (1 - p) ./ (prec + 1)) * scaler;
  }
}
