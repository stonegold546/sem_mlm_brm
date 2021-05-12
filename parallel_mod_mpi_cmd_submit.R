data.location <- "summary_text/mpi_dat/"  # location of csv dataset

# Load packages (Check for packages that should be installed) ----
source("parallel_mod_mpi_cmd_load_packages.R")

# Pre-process data (Worth reviewing source file) ----
source("parallel_mod_mpi_cmd_prep_data.R")

stan.scripts <- "summary_text/stan_scripts/"  # location of stan scripts
print.images <- "summary_text/text_clean/"  # location to print images

# FIGURE 1 ----
ggplot(X, aes(response.l)) + geom_histogram(col = 1) + facet_wrap(~ item_n) +
  scale_x_continuous(name = "Logit-prevalence of deprivation") +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(.5, "cm"), panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(y = "Number of countries")
ggsave(paste0(print.images, "01_deprivation.pdf"), width = 6, height = 2.75)

# MODELLING ----

# CONGENERIC CFA ----
# Congeneric CFA in lavaan ----
# Creating error SD as square root of error variance
writeLines(cong.lav.form <- paste(
  paste0("F =~ ", paste0(colnames(dat)[15:20], collapse = " + ")),
  "x1 ~~ a * x1\nx2 ~~ b * x2\nx3 ~~ c * x3\nx4 ~~ d * x4\nx5 ~~ e * x5\nx6 ~~ f * x6",
  "a1 := sqrt(a)\nb1 := sqrt(b)\nc1 := sqrt(c)\nd1 := sqrt(d)\ne1 := sqrt(e)\nf1 := sqrt(f)",
  sep = "\n"))
summary(cong.cfa <- cfa(
  cong.lav.form, dat, std.lv = TRUE, meanstructure = TRUE, missing = "ML"),
  fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
parameterEstimates(cong.cfa, standardized = TRUE, rsquare = TRUE)
#    lhs op     rhs label    est    se       z pvalue ci.lower ci.upper std.lv std.all std.nox
# 1    F =~      x1        2.389 0.215  11.088  0.000    1.967    2.811  2.389   0.875   0.875
# 2    F =~      x2        1.650 0.139  11.838  0.000    1.376    1.923  1.650   0.908   0.908
# 3    F =~      x3        1.137 0.120   9.492  0.000    0.903    1.372  1.137   0.790   0.790
# 4    F =~      x4        2.384 0.177  13.457  0.000    2.036    2.731  2.384   0.975   0.975
# 5    F =~      x5        1.573 0.149  10.550  0.000    1.281    1.866  1.573   0.848   0.848
# 6    F =~      x6        1.884 0.156  12.063  0.000    1.578    2.190  1.884   0.918   0.918
# 14  x1 ~1                0.053 0.272   0.194  0.846   -0.480    0.586  0.053   0.019   0.019
# 15  x2 ~1               -0.579 0.181  -3.204  0.001   -0.933   -0.225 -0.579  -0.319  -0.319
# 16  x3 ~1               -1.425 0.143  -9.948  0.000   -1.706   -1.144 -1.425  -0.990  -0.990
# 17  x4 ~1               -1.667 0.243  -6.849  0.000   -2.145   -1.190 -1.667  -0.682  -0.682
# 18  x5 ~1               -0.696 0.185  -3.764  0.000   -1.058   -0.334 -0.696  -0.375  -0.375
# 19  x6 ~1               -2.498 0.204 -12.234  0.000   -2.899   -2.098 -2.498  -1.217  -1.217
# 21  a1 := sqrt(a)    a1  1.321 0.103  12.843  0.000    1.120    1.523  1.321   0.484   0.484
# 22  b1 := sqrt(b)    b1  0.760 0.063  12.164  0.000    0.638    0.883  0.760   0.419   0.419
# 23  c1 := sqrt(c)    c1  0.882 0.065  13.509  0.000    0.754    1.010  0.882   0.613   0.613
# 24  d1 := sqrt(d)    d1  0.546 0.084   6.472  0.000    0.381    0.712  0.546   0.223   0.223
# 25  e1 := sqrt(e)    e1  0.983 0.075  13.117  0.000    0.836    1.130  0.983   0.530   0.530
# 26  f1 := sqrt(f)    f1  0.814 0.068  11.892  0.000    0.680    0.948  0.814   0.397   0.397
# 27  x1 r2      x1        0.766    NA      NA     NA       NA       NA     NA      NA      NA
# 28  x2 r2      x2        0.825    NA      NA     NA       NA       NA     NA      NA      NA
# 29  x3 r2      x3        0.624    NA      NA     NA       NA       NA     NA      NA      NA
# 30  x4 r2      x4        0.950    NA      NA     NA       NA       NA     NA      NA      NA
# 31  x5 r2      x5        0.719    NA      NA     NA       NA       NA     NA      NA      NA
# 32  x6 r2      x6        0.843    NA      NA     NA       NA       NA     NA      NA      NA

# Congeneric regression in Stan ----
# Data list with references to equation 11
data <- list(
  alpha_scale = 5 / qnorm(.975),  # scale of a parameter
  beta_scale = 2 / qnorm(.975),  # scale of sigma_t
  sigma_scale = 5 / qnorm(.975),  # scale of sigma_i
  lambda_median = log(.5),  # log-median of lambda_i
  lambda_scale = log(4 / .5) / qnorm(.975),  # scale of lambda_i
  # item id number long-form, ids be contiguous from 1 to number of items
  item_id = X$item,
  Ni = max(X$item),  # number of items
  # respondent id number long-form, ids be contiguous from 1 to number of respondents
  resp_id = X$id,
  Np = max(X$id),  # number of respondents
  y = X$response.l,  # logit-response
  N = length(X$response),  # rows in data (long form)
  ret_yhat = 1,  # return predicted outcomes (1: yes, 0: no)
  ret_ll = 0  # return log-likelihood (1: yes, 0 = no),
)

X.tmp <- lavaan::HolzingerSwineford1939[, paste0("x", 1:3)]
colnames(X.tmp) <- paste0("x.", 1:ncol(X.tmp))
head(X.tmp <- reshape(as.data.frame(scale(X.tmp, center = FALSE)), direction = "long", 1:ncol(X.tmp)))

data <- list(
  alpha_scale = 5 / qnorm(.975),  # scale of a parameter
  beta_scale = 2 / qnorm(.975),  # scale of sigma_t
  sigma_scale = 5 / qnorm(.975),  # scale of sigma_i
  lambda_median = log(.5),  # log-median of lambda_i
  lambda_scale = log(4 / .5) / qnorm(.975),  # scale of lambda_i
  # item id number long-form, ids be contiguous from 1 to number of items
  item_id = X.tmp$time,
  Ni = max(X.tmp$time),  # number of items
  # respondent id number long-form, ids be contiguous from 1 to number of respondents
  resp_id = X.tmp$id,
  Np = max(X.tmp$id),  # number of respondents
  y = X.tmp$x,  # logit-response
  N = length(X.tmp$x),  # rows in data (long form)
  ret_yhat = 1,  # return predicted outcomes (1: yes, 0: no)
  ret_ll = 0  # return log-likelihood (1: yes, 0 = no),
)

# Compile Stan model (takes time)
cong.mod <- cmdstan_model(file.path(paste0(stan.scripts, "cong.stan")))
cong.fit <- cong.mod$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems,
)
cong.fit$cmdstan_diagnose()  # Run quick diagnostic checks
cong.fit.rs <- read_stan_csv(cong.fit$output_files())  # convert to Rstan

# Print major parameters
print(cong.fit.rs, c("theta_p", "yhat", "lp__"), include = FALSE)
#             mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha      -1.10    0.01 0.46 -2.01 -1.39 -1.11 -0.82 -0.19  4576 1.00
# sigma_beta  1.02    0.00 0.32  0.57  0.79  0.96  1.18  1.82 15613 1.00
# beta[1]     0.03    0.01 0.26 -0.48 -0.15  0.03  0.20  0.54   900 1.01
# beta[2]    -0.58    0.01 0.17 -0.92 -0.70 -0.58 -0.47 -0.24   817 1.01
# beta[3]    -1.42    0.00 0.14 -1.69 -1.52 -1.42 -1.33 -1.15  1078 1.01
# beta[4]    -1.67    0.01 0.23 -2.11 -1.82 -1.67 -1.51 -1.21   716 1.01
# beta[5]    -0.70    0.01 0.18 -1.05 -0.82 -0.70 -0.58 -0.35   965 1.01
# beta[6]    -2.49    0.01 0.20 -2.87 -2.62 -2.49 -2.36 -2.10   815 1.01
# lambda[1]   2.33    0.00 0.21  1.94  2.18  2.32  2.47  2.77  1927 1.00
# lambda[2]   1.61    0.00 0.14  1.37  1.52  1.61  1.70  1.90  1643 1.00
# lambda[3]   1.11    0.00 0.12  0.89  1.03  1.10  1.18  1.35  2619 1.00
# lambda[4]   2.34    0.00 0.17  2.03  2.22  2.33  2.45  2.71  1316 1.00
# lambda[5]   1.53    0.00 0.15  1.26  1.43  1.53  1.63  1.84  2046 1.00
# lambda[6]   1.84    0.00 0.15  1.56  1.74  1.84  1.94  2.16  1618 1.00
# sigma[1]    1.35    0.00 0.11  1.16  1.28  1.35  1.42  1.59 19634 1.00
# sigma[2]    0.78    0.00 0.07  0.66  0.73  0.78  0.82  0.92 14844 1.00
# sigma[3]    0.90    0.00 0.07  0.78  0.86  0.90  0.95  1.05 19643 1.00
# sigma[4]    0.55    0.00 0.09  0.37  0.50  0.55  0.61  0.73  3596 1.00
# sigma[5]    1.01    0.00 0.08  0.87  0.95  1.00  1.06  1.18 19595 1.00
# sigma[6]    0.84    0.00 0.07  0.70  0.79  0.83  0.88  0.99 13708 1.00

# saveRDS(cong.fit.rs, "res/mod_1.rds")
cong.fit.rs <- readRDS("res/mod_1.rds")

# Traceplots:
mcmc_trace(cong.fit.rs, regex_pars = c("alpha", "beta", "lambda", "sigma"))
# Rank plots:
mcmc_rank_overlay(
  cong.fit.rs, regex_pars = c("alpha", "beta", "lambda", "sigma"),
  ref_line = TRUE)

cong.pars.df <- as.data.frame(
  summary(cong.fit.rs, c("beta", "lambda", "sigma"))$summary)
# beta[1-6] are means
# lambda[1-6] are loadings
# sigma[1-6] are residual standard deviations
cong.pars.df$variable <- rownames(cong.pars.df)
# Join lavaan and Bayes results
cong.pars.df
cong.pars.df.joined <- as.data.table(cbind(
  cong.pars.df[, c(11, 6, 4, 8)],
  parameterEstimates(cong.cfa)[c(14:19, 1:6, 21:26), c(5, 9:10)]))
cong.pars.df.joined
# Rename columns
setnames(cong.pars.df.joined, 2:7, c("m.1", "ll.1", "ul.1", "m.2", "ll.2", "ul.2"))
cong.pars.df.joined <- reshape(
  cong.pars.df.joined, direction = "long", varying = list(c(2, 5), c(3, 6), c(4, 7)))
cong.pars.df.joined[, model := c("Bayesian", "ML")[time]]
cong.pars.df.joined[, type := factor(
  gsub("\\[\\d\\]", "", variable), levels = c("lambda", "beta", "sigma"),
  labels = c("Loadings", "Means", "Residual SD"))]
cong.pars.df.joined[, item := as.integer(gsub("[a-z]|\\[|\\]", "", variable))]
cong.pars.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
ggplot(cong.pars.df.joined, aes(item, m.1, fill = model, shape = model)) +
  geom_pointrange(aes(ymin = ll.1, ymax = ul.1), position = position_dodge(.5), alpha = .65) +
  theme_classic() + theme(legend.position = "top", strip.background = element_blank(),
                          axis.title.y = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ type, scales = "free_y") +
  labs(x = "Indicator number", subtitle = "Estimates with 95% interval",
       shape = "Method", fill = "Method", caption = fig.capt)
ggsave(paste0(print.images, "02_cong_params_lav_bayes.pdf"), width = 6.5, height = 4)

# Extract predicted outcomes: 500 samples randomly drawn from 24000 samples
dim(G.cong.f <- as.data.frame(cong.fit.rs, "yhat"))
# Should be 100 samples by 603 responses
set.seed(12345)
dim(G.cong <- G.cong.f[sample(1:nrow(G.cong.f), 100), ])
G.cong <- data.table(item = X$item_n, response = X$response.l, t(G.cong))
setnames(G.cong, 3:ncol(G.cong), paste0("post.", 1:(ncol(G.cong) - 2)))
G.cong <- reshape(G.cong, direction = "long", varying = 3:ncol(G.cong))
G.cong

ggplot(G.cong, aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "#999999") +
  geom_density(col = 1, fill = NA, data = G.cong[time == 1]) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

# Compile Stan model (takes time)
cong.mod.j <- cmdstan_model(file.path(paste0(stan.scripts, "cong_joint.stan")))
cong.fit.j <- cong.mod.j$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems
)
cong.fit.j$cmdstan_diagnose()  # Run quick diagnostic checks
cong.fit.j.rs <- read_stan_csv(cong.fit.j$output_files())  # convert to Rstan

# Print major parameters
print(cong.fit.j.rs, c("Theta", "yhat", "lp__"), include = FALSE)
#                mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha         -1.10    0.01 0.46 -2.00 -1.39 -1.11 -0.82 -0.17  6938    1
# sigma_beta     0.99    0.00 0.31  0.55  0.77  0.93  1.15  1.76 19410    1
# beta[1]       -0.01    0.01 0.27 -0.54 -0.19 -0.01  0.17  0.51  1634    1
# beta[2]       -0.58    0.00 0.18 -0.93 -0.70 -0.58 -0.46 -0.24  1598    1
# beta[3]       -1.40    0.00 0.14 -1.69 -1.50 -1.40 -1.31 -1.13  1914    1
# beta[4]       -1.68    0.01 0.24 -2.15 -1.84 -1.68 -1.52 -1.22  1383    1
# beta[5]       -0.71    0.00 0.18 -1.07 -0.83 -0.71 -0.59 -0.36  1807    1
# beta[6]       -2.44    0.00 0.19 -2.81 -2.57 -2.44 -2.31 -2.08  1701    1
# lambda[1]      2.40    0.00 0.21  2.02  2.26  2.39  2.54  2.84  3140    1
# lambda[2]      1.59    0.00 0.13  1.35  1.50  1.59  1.68  1.87  2890    1
# lambda[3]      1.17    0.00 0.11  0.96  1.09  1.17  1.24  1.41  3589    1
# lambda[4]      2.34    0.00 0.17  2.03  2.22  2.33  2.45  2.70  2263    1
# lambda[5]      1.53    0.00 0.15  1.27  1.43  1.53  1.63  1.83  3564    1
# lambda[6]      1.69    0.00 0.15  1.42  1.58  1.68  1.79  2.00  2974    1
# ln_alpha      -0.42    0.00 0.36 -1.15 -0.61 -0.41 -0.22  0.29 10051    1
# ln_sigma_beta  0.76    0.00 0.39  0.32  0.51  0.67  0.90  1.77  9871    1
# ln_sigma2[1]   0.36    0.00 0.18  0.00  0.24  0.37  0.49  0.72 10894    1
# ln_sigma2[2]  -0.51    0.00 0.16 -0.83 -0.62 -0.51 -0.40 -0.19 21408    1
# ln_sigma2[3]  -0.74    0.00 0.21 -1.17 -0.88 -0.74 -0.60 -0.34  5860    1
# ln_sigma2[4]  -1.11    0.00 0.29 -1.72 -1.29 -1.10 -0.92 -0.60  6632    1
# ln_sigma2[5]  -0.01    0.00 0.15 -0.30 -0.12 -0.02  0.09  0.29 23470    1
# ln_sigma2[6]  -0.51    0.00 0.19 -0.89 -0.63 -0.50 -0.38 -0.14  9014    1
# ln_lambda[1]   0.57    0.00 0.20  0.20  0.43  0.56  0.69  0.99 11663    1
# ln_lambda[2]   0.16    0.00 0.10  0.03  0.09  0.14  0.21  0.41 20368    1
# ln_lambda[3]   1.07    0.00 0.23  0.66  0.91  1.05  1.21  1.54  9260    1
# ln_lambda[4]   0.26    0.00 0.16  0.05  0.14  0.23  0.35  0.67 13097    1
# ln_lambda[5]   0.15    0.00 0.09  0.03  0.08  0.13  0.20  0.38 20001    1
# ln_lambda[6]   0.62    0.00 0.25  0.18  0.44  0.61  0.78  1.14  7527    1
# r_tr           0.12    0.00 0.05  0.03  0.08  0.11  0.15  0.23  1767    1
# r             -0.77    0.00 0.10 -0.94 -0.85 -0.78 -0.70 -0.53  1767    1

# saveRDS(cong.fit.j.rs, "res/mod_2.rds")
cong.fit.j.rs <- readRDS("res/mod_2.rds")

# Traceplots:
mcmc_trace(cong.fit.j.rs, regex_pars = c("beta", "lambda", "sigma2"))
# Rank plots:
mcmc_rank_overlay(
  cong.fit.j.rs, regex_pars = c("beta", "lambda", "sigma2"),
  ref_line = TRUE)

G.cong.f.j <- as.data.frame(cong.fit.j.rs, "yhat")
# Should be 100 samples by 603 responses
set.seed(12345)
dim(G.cong.j <- G.cong.f.j[sample(1:nrow(G.cong.f.j), 100), ])
G.cong.j <- data.table(item = X$item_n, response = X$response.l, t(G.cong.j))
setnames(G.cong.j, 3:ncol(G.cong.j), paste0("post.", 1:(ncol(G.cong.j) - 2)))
G.cong.j <- reshape(G.cong.j, direction = "long", varying = 3:ncol(G.cong.j))
G.cong.j

ggplot(G.cong.j, aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "darkgrey") +
  geom_density(aes(response), col = 1, fill = NA, data = G.cong.j[time == 1]) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

G.cong.new <- rbindlist(list(G.cong, G.cong.j), idcol = "model")
G.cong.new[, model.t := c("Location", "Location-scale")[model]]

ggplot(G.cong.new, aes(response)) +
  facet_grid_paginate(item ~ model.t, scales = "free_y") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "darkgrey", size = .25) +
  geom_density(col = 1, fill = NA, data = G.cong.j[time == 1]) +
  coord_cartesian(xlim = c(-10, 10)) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))
ggsave(paste0(print.images, "03_cong_ppc.pdf"), width = 6.5, height = 6.5)
rm(G.cong.new)

G.cong.new <- data.table(
  m1 = sapply(1:ncol(G.cong.f), function (i) median(abs(G.cong.f[, i] - X$response.l[i]))),
  m2 = sapply(1:ncol(G.cong.f.j), function (i) median(abs(G.cong.f.j[, i] - X$response.l[i])))
)

ggplot(G.cong.new, aes(m1, m2)) + coord_fixed() +
  geom_point(shape = 4, alpha = .25, size = .25) +
  geom_abline(col = "#E69F00", alpha = 1) +
  geom_smooth(formula = y ~ s(x, bs = "ts"), method = "gam", size = .5) +
  theme_classic() +
  labs(x = "Location error", y = "Location-scale error",
       caption = "Straight line is 45-degree line. Smooth fit is GAM fit")
ggsave(paste0(print.images, "04_compare_mad.pdf"), width = 4, height = 4)

describe(G.cong.new)
#    vars   n mean   sd median trimmed  mad  min  max range skew kurtosis   se
# m1    1 603 0.90 0.46   0.73    0.81 0.25 0.45 3.99  3.54 2.27     6.88 0.02
# m2    2 603 0.86 0.49   0.72    0.77 0.29 0.24 4.10  3.87 2.16     6.33 0.02

cong.pars.j.df.2 <- as.data.frame(summary(
  cong.fit.j.rs, c("beta", "lambda", "ln_sigma2", "ln_lambda", "Theta"))$summary)
cong.pars.j.df.2 <- cong.pars.j.df.2[, c(1, 6)]
cong.pars.j.df.2$variable <- rownames(cong.pars.j.df.2)
cong.pars.j.df.2 <- as.data.table(cong.pars.j.df.2)
# Use equation to create predicted outcome in wide form
cong.pars.j.df.2.res <-
  t(sapply(cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",1]"), mean], function (theta) {
    theta * cong.pars.j.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean] +
      cong.pars.j.df.2[variable %in% paste0("beta[", 1:6, "]"), mean]
  }))
colnames(cong.pars.j.df.2.res) <- paste0("x.", 1:6)
# Add standardized averages and row id
cong.pars.j.df.2.res <- data.frame(
  cong.pars.j.df.2.res,
  theta = cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",1]"), mean], id = 1:101)
# Transform data to long form
cong.pars.j.df.2.res <- as.data.table(
  reshape(cong.pars.j.df.2.res, direction = "long", varying = 1:6, timevar = "item"))
# Add lambda parameter
cong.pars.j.df.2.res[, lambda := cong.pars.j.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean][item]]
plt.1 <- ggplot(cong.pars.j.df.2.res, aes(theta, x, group = item, label = item)) +
  geom_line(aes(alpha = lambda)) + geom_rug(y = NA, alpha = .1) +
  geom_dl(aes(alpha = lambda), method = "last.points",
          data = cong.pars.j.df.2.res[item %in% c(1, 2, 3, 5)]) +
  geom_dl(aes(alpha = lambda), method = "first.points",
          data = cong.pars.j.df.2.res[item %in% c(1, 4, 5, 6)]) +
  scale_alpha(range = c(.3, 1)) + scale_size(range = c(.25, 1)) + guides(alpha = FALSE, size = FALSE) +
  # scale_y_continuous(labels = percent_format(), breaks = seq(0, 4, .2)) +
  theme(panel.grid.minor.y = element_blank()) +
  labs(y = "Expected average", x = "Country location")
plt.1

cong.pars.j.df.2.ln.res <-
  t(sqrt(exp(sapply(cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",2]"), mean], function (theta) {
    theta * cong.pars.j.df.2[variable %in% paste0("ln_lambda[", 1:6, "]"), mean] +
      cong.pars.j.df.2[variable %in% paste0("ln_sigma2[", 1:6, "]"), mean]
  }))))
colnames(cong.pars.j.df.2.ln.res) <- paste0("x.", 1:6)
# Add standardized averages and row id
cong.pars.j.df.2.ln.res <- data.frame(
  cong.pars.j.df.2.ln.res,
  theta = cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",2]"), mean], id = 1:101)
# Transform data to long form
cong.pars.j.df.2.ln.res <- as.data.table(
  reshape(cong.pars.j.df.2.ln.res, direction = "long", varying = 1:6, timevar = "item"))
# Add lambda parameter
cong.pars.j.df.2.ln.res[, lambda := cong.pars.j.df.2[variable %in% paste0("ln_lambda[", 1:6, "]"), mean][item]]
plt.2 <- ggplot(cong.pars.j.df.2.ln.res, aes(theta, x, group = item, label = item)) +
  geom_line(aes(alpha = lambda)) + geom_rug(y = NA, alpha = .1) +
  geom_dl(aes(alpha = lambda), method = "last.points") +
  scale_alpha(range = c(.3, 1)) + scale_size(range = c(.25, 1)) + guides(alpha = FALSE, size = FALSE) +
  # scale_y_continuous(labels = percent_format(), breaks = seq(0, 4, .2)) +
  theme(panel.grid.minor.y = element_blank()) +
  labs(y = "Expected SD", x = "Country log-variance",
       caption = fig.capt)
plt.2

plt.1 + plt.2

cong.pars.j.df.3.res <- data.table(
  loc = cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",1]"), mean],
  scale = cong.pars.j.df.2[variable %in% paste0("Theta[", 1:101, ",2]"), mean])
plt.3 <- ggplot(cong.pars.j.df.3.res, aes(loc, scale)) + geom_point(shape = 4) +
  labs(x = "Country location", y = "Country log-variance")
plt.3

(plt.1 + plt.2) / plt.3 + plot_layout(heights = c(.7, .3))
ggsave(paste0(print.images, "05_model_suggestions_1.pdf"), width = 6.5, height = 6)

colnames(cong.pars.j.df.2.ln.res)[c(1, 4:5)] <- paste0("ln_", colnames(cong.pars.j.df.2.ln.res)[c(1, 4:5)])
cong.pars.res.df <- merge(cong.pars.j.df.2.res, cong.pars.j.df.2.ln.res)
cong.pars.res.df <- merge(
  cong.pars.res.df, aggregate(id ~ item + item_n, X, length)[, -3], all.x = TRUE, by = "item")
cong.pars.res.df

ggplot(cong.pars.res.df, aes(theta, x)) +
  # geom_line(aes(y = qnorm(.025, x, sqrt(ln_x))), alpha = .25) +
  # geom_line(aes(y = qnorm(.975, x, sqrt(ln_x))), alpha = .25) +
  geom_ribbon(aes(ymin = qnorm(.025, x, sqrt(ln_x)),
                  ymax = qnorm(.975, x, sqrt(ln_x))), alpha = .25) +
  geom_line() +
  facet_wrap(~ item_n, scales = "free_y") + theme_classic() +
  theme(strip.background = element_blank()) +
  labs(x = "Country location", y = "Country prediction with 95% interval")
ggsave(paste0(print.images, "06_model_suggestions_2.pdf"), width = 6.5, height = 4)

# On probability scale (unpresented)
ggplot(cong.pars.res.df, aes(theta, plogis(x))) + geom_line() +
  geom_ribbon(aes(ymin = plogis(qnorm(.025, x, sqrt(ln_x))),
                  ymax = plogis(qnorm(.975, x, sqrt(ln_x)))), alpha = .25) +
  facet_wrap(~ item) + theme_classic()

(ppc <- ppc.j <- data.table(item = 1:6))

est.function <- mean
(est.s <- aggregate(response.l ~ item, X, est.function)[, 2])
# [1]  0.08278299 -0.57900871 -1.42493528 -1.66909220 -0.68003053 -2.49829869
percent(ppc$mean <- sapply(1:6, function (i) {
  mean(apply(G.cong.f[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "-5.33%" "-1.25%" "1.15%"  "0.37%"  "-1.80%" "4.57%" 
percent(ppc.j$mean <- sapply(1:6, function (i) {
  mean(apply(G.cong.f.j[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "-11.55%" "1.80%"   "7.59%"   "-1.77%"  "-2.76%"  "19.07%" 

est.function <- median
(est.s <- aggregate(response.l ~ item, X, est.function)[, 2])
# [1]  0.2735150 -0.3380821 -1.0833937 -1.4121223 -0.5277053 -1.9025473
percent(ppc$median <- sapply(1:6, function (i) {
  mean(apply(G.cong.f[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "-12.3%" "-28.7%" "-46.7%" "-4.0%"  "-16.7%" "-49.0%"
percent(ppc.j$median <- sapply(1:6, function (i) {
  mean(apply(G.cong.f.j[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "5.4%"   "-24.5%" "-21.4%" "-4.1%"  "-11.5%" "-46.4%"

est.function <- sd
(est.s <- aggregate(response.l ~ item, X, est.function)[, 2])
# [1] 2.737565 1.825428 1.446722 2.468665 1.864159 2.062517
percent(ppc$sd <- sapply(1:6, function (i) {
  mean(apply(G.cong.f[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "-2.633%" "-1.000%" "-0.871%" "-0.812%" "-1.817%" "-1.442%"
percent(ppc.j$sd <- sapply(1:6, function (i) {
  mean(apply(G.cong.f.j[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "9.3%"   "-4.9%"  "17.2%"  "8.1%"   "-0.1%"  "-33.8%"

est.function <- mad
(est.s <- aggregate(response.l ~ item, X, est.function)[, 2])
# [1] 2.761740 1.950652 1.314030 3.310832 2.044396 2.062330
percent(ppc$mad <- sapply(1:6, function (i) {
  mean(apply(G.cong.f[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "31.5%"  "26.9%"  "40.0%"  "-41.5%" "-0.8%"  "43.5%" 
percent(ppc.j$mad <- sapply(1:6, function (i) {
  mean(apply(G.cong.f.j[, which(X$item == i)], 1, est.function) > est.s[i])
}) - .5)
# [1] "32.4%"  "22.4%"  "35.4%"  "-42.4%" "-1.3%"  "14.8%" 

(ppc.res <- rbindlist(list(ppc, ppc.j), idcol = "model"))
ppc.res[, model.t := c("Basic", "Location-Scale")[model]]
ppc.res

setnames(ppc.res, 3:6, paste0("value.", 1:4))

(ppc.res <- reshape(ppc.res, direction = "long", varying = 3:6, timevar = "metric"))
ppc.res[, metric.t := c("mean", "median", "sd", "mad")[metric]]
ppc.res[, loc := c("scale", "location")[(metric <= 2) + 1]]
ppc.res[, rob := c("non-robust", "robust")[(metric %% 2 == 0) + 1]]
ppc.res

ggplot(ppc.res, aes(value + .0, reorder(item, -item))) +
  geom_point(shape = c("L", "S")[ppc.res$model]) +
  geom_vline(xintercept = .0, linetype = 2, alpha = .25) +
  facet_wrap(~ reorder(metric.t, metric)) +
  scale_x_continuous(labels = percent_format()) +
  theme(legend.position = "top")

ppc.dist <- data.table(x = X$response.l)
ppc.dist <- cbind(ppc.dist, t(apply(G.cong, 2, quantile, c(.25, .75))))
ppc.dist <- cbind(ppc.dist, t(apply(G.cong.j, 2, quantile, c(.25, .75))))
setnames(ppc.dist, 2:5, c("ll", "ul", "llj", "ulj"))
ppc.dist[order(x), id := 1:.N]
ppc.dist

ppc.dist[, mean(x > ll & x < ul)]
ppc.dist[, mean(x > llj & x < ulj)]

ggplot(ppc.dist, aes(id, x)) +
  geom_linerange(aes(ymin = ll, ymax = ul), alpha = .25) +
  geom_point()
ggplot(ppc.dist, aes(id, x)) +
  geom_linerange(aes(ymin = llj, ymax = ulj), alpha = .25) +
  geom_point()

# R Session Info (Contains exact package versions and system information) ----
sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 20.10
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] ggforce_0.3.2          cmdstanr_0.1.0         rstan_2.21.2           StanHeaders_2.21.0-5   bayesplot_1.7.2       
#  [6] ggrepel_0.8.2          directlabels_2020.6.17 ggplot2_3.3.2          latex2exp_0.4.0        patchwork_1.0.1       
# [11] scales_1.1.1           tidyr_1.1.2            dplyr_1.0.2            glmmTMB_1.0.2.1        lavaan_0.6-6          
# [16] psych_2.0.7            data.table_1.13.0     
# 
# loaded via a namespace (and not attached):
#  [1] nlme_3.1-149       matrixStats_0.56.0 bit64_0.9-7.1      backports_1.1.10   tools_4.0.3        TMB_1.7.18        
#  [7] utf8_1.1.4         R6_2.4.1           colorspace_1.4-1   withr_2.3.0        tidyselect_1.1.0   gridExtra_2.3     
# [13] prettyunits_1.1.1  mnormt_2.0.1       processx_3.4.4     emmeans_1.5.1      bit_1.1-15.2       curl_4.3          
# [19] compiler_4.0.3     cli_2.1.0          sandwich_2.5-1     posterior_0.1.2    labeling_0.4.2     checkmate_2.0.0   
# [25] mvtnorm_1.1-1      quadprog_1.5-8     ggridges_0.5.2     callr_3.5.1        stringr_1.4.0      digest_0.6.27     
# [31] pbivnorm_0.6.0     minqa_1.2.4        pkgconfig_2.0.3    lme4_1.1-23        rlang_0.4.8        rstudioapi_0.11   
# [37] generics_0.0.2     farver_2.0.3       zoo_1.8-8          jsonlite_1.7.1     vroom_1.2.1        inline_0.3.15     
# [43] magrittr_1.5       loo_2.3.1          Matrix_1.2-18      Rcpp_1.0.5         munsell_0.5.0      fansi_0.4.1       
# [49] abind_1.4-5        lifecycle_0.2.0    stringi_1.5.3      multcomp_1.4-13    MASS_7.3-53        pkgbuild_1.1.0    
# [55] plyr_1.8.6         grid_4.0.3         parallel_4.0.3     crayon_1.3.4       lattice_0.20-41    splines_4.0.3     
# [61] tmvnsim_1.0-2      ps_1.4.0           pillar_1.4.6       boot_1.3-25        estimability_1.3   reshape2_1.4.4    
# [67] codetools_0.2-16   stats4_4.0.3       glue_1.4.2         V8_3.3.0           RcppParallel_5.0.2 tweenr_1.0.1      
# [73] vctrs_0.3.4        nloptr_1.2.2.2     polyclip_1.10-0    gtable_0.3.0       purrr_0.3.4        assertthat_0.2.1  
# [79] xtable_1.8-4       coda_0.19-3        survival_3.2-7     tibble_3.0.4       statmod_1.4.34     TH.data_1.0-10    
# [85] ellipsis_0.3.1    
