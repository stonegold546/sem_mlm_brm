// congeneric with log-normal loadings
// Distribution is beta, mean & sample size parameterization
// Cov-mat decomp approach
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
  vector<lower = 0>[Ni] prec;
  vector[Ni] latent_mv[Np];
}
model {
  alpha ~ normal(0, alpha_scale);
  sigma_beta ~ normal(0, beta_scale);
  beta ~ normal(alpha, sigma_beta);

  lambda ~ lognormal(lambda_median, lambda_scale);

  prec ~ gamma(2, .1);

  {
    matrix[Ni, Ni] Sigma = tcrossprod(to_matrix(lambda)) + diag_matrix(rep_vector(.01, Ni));
    vector[N] prob;

    latent_mv ~ multi_normal(beta, Sigma);
    for (i in 1:N) prob[i] = inv_logit(latent_mv[resp_id[i], item_id[i]]);
    y ~ beta_proportion(prob, prec[item_id]);
  }
}
generated quantities {
  vector[Ni] i_means;
  vector[Ni] sigma;

  {
    vector[Ni] p = inv_logit(beta);
    i_means = p * scaler;
    sigma = sqrt(p .* (1 - p) ./ (prec + 1)) * scaler;
  }
}
