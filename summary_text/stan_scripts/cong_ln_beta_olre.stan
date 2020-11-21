// congeneric with log-normal loadings & OLRE
// Distribution is beta, mean & sample size parameterization
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
  real<lower = 0> prec;
  real<lower = 0> nu;
  real<lower = 0> tau;
  vector<lower = 0>[Ni] olre_sd;
  vector[N] olre;
}
model {
  alpha ~ normal(0, alpha_scale);
  sigma_beta ~ normal(0, beta_scale);
  beta ~ normal(alpha, sigma_beta);

  lambda ~ lognormal(lambda_median, lambda_scale);
  theta_p ~ std_normal();

  olre_sd ~ student_t(3, 0, 1);
  nu ~ gamma(2, .1);
  tau ~ gamma(nu / 2, nu / 2);
  olre ~ std_normal();

  prec ~ gamma(2, .1);

  {
    vector[N] prob = inv_logit(
      beta[item_id] + lambda[item_id] .* theta_p[resp_id] +
      olre_sd[item_id] .* olre / sqrt(tau));
    y ~ beta_proportion(prob, prec);
  }
}
generated quantities {
  vector[Ni] i_means = inv_logit(beta) * scaler;
  vector[Nll] log_lik;
  vector[Ny] yhat;

  {
    real p;
    for (i in 1:max(Nll, Ny)) {
      p = inv_logit(
        beta[item_id[i]] + lambda[item_id[i]] * theta_p[resp_id[i]] +
        olre_sd[item_id[i]] * olre[i] / sqrt(tau));
      if (Nll > 0) log_lik[i] = beta_proportion_lpdf(y[i] | p, prec);
      if (Ny > 0) yhat[i] = beta_proportion_rng(p, prec) * scaler;
    }
  }
}
