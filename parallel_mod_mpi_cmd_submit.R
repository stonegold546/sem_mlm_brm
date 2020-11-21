data.location <- "summary_text/mpi_dat/"  # location of csv dataset

# Load packages (Check for packages that should be installed) ----
source("parallel_mod_mpi_cmd_load_packages.R")

# Pre-process data (Worth reviewing source file) ----
source("parallel_mod_mpi_cmd_prep_data.R")

stan.scripts <- "summary_text/stan_scripts/"  # location of stan scripts
print.images <- "summary_text/text_clean/"  # location to print images

# FIGURE 1 ----
X %>%
  ggplot(aes(response)) + geom_histogram(col = 1) + facet_wrap(~ item_n) +
  scale_x_continuous(labels = percent_format(), name = "Prevalence of deprivation (%)") +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(.5, "cm"), panel.grid.minor = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = "Number of countries")
ggsave(paste0(print.images, "01_deprivation.pdf"), width = 6, height = 3)

# MODELLING ----

# Parallel CFA ----
# Parallel CFA in lavaan ----
writeLines(para.lav.form <- paste(
  paste("F =~ a *", paste0(colnames(dat)[15:20], collapse = " + a * ")),
  paste0("x", 1:6, " ~~ b * x", 1:6, collapse = "\n"), sep = "\n"
))
summary(para.cfa <- cfa(
  para.lav.form, dat, std.lv = TRUE, meanstructure = TRUE, missing = "ML"
), fit.measures = TRUE)

# Parallel CFA as regression ----
(para.reg <- glmmTMB(response ~ 0 + item_t + (1 | id), X))

# Compare (parameter) estimates from lavaan and glmmTMB ----
logLik(para.cfa); logLik(para.reg)  # log-lieklihoods
# Means
as.matrix(parameterEstimates(para.cfa)[14:19, 5:7])
coef(summary(para.reg))$cond
# Constrained error variance
as.matrix(parameterEstimates(para.cfa)[7, 5:7])
sigma(para.reg) ^ 2
# Constrained loading
as.matrix(parameterEstimates(para.cfa)[1, 5:7])
attr(VarCorr(para.reg)$cond$id, "stddev")

# Tau-equivalent CFA ----
# Tau-equivalent CFA in lavaan ----
writeLines(te.lav.form <- paste(
  "F =~ a *", paste0(colnames(dat)[15:20], collapse = " + a * ")
))
summary(te.cfa <- cfa(
  te.lav.form, dat, std.lv = TRUE, meanstructure = TRUE, missing = "ML"
), fit.measures = TRUE)

# Tau-equivalent CFA as regression ----
(te.reg <- glmmTMB(
  response ~ 0 + item_t + (1 | id), X,
  # Below add indicators as predictors of variance
  dispformula = ~ 0 + item_t))

# Compare (parameter) estimates from lavaan and glmmTMB ----
logLik(te.cfa); logLik(te.reg)  # log-lieklihoods
# Means
as.matrix(parameterEstimates(te.cfa)[14:19, 5:7])
coef(summary(te.reg))$cond
# Error variances
as.matrix(parameterEstimates(te.cfa)[7:12, 5:7])
exp(fixef(te.reg)$disp)
# Constrained loading
as.matrix(parameterEstimates(te.cfa)[1, 5:7])
attr(VarCorr(te.reg)$cond$id, "stddev")

# Model comparing parallel to tau-equivalent ----
anova(para.reg, te.reg)
anova(para.cfa, te.cfa)

# FIGURE 2 ----
X$te.preds <- predict(te.reg)  # Save predictions
# Extract random effects and standardize them
X$ave.te <- (ranef(te.reg)$cond$id[, 1])[X$id] / attr(VarCorr(te.reg)$cond$id, "stddev")
X$item_abb <- paste0(substr(X$item_n, 1, 5), ".")  # Create abbreviations for item names

X %>%
  ggplot(aes(ave.te, fitted(te.reg), group = item_t, label = item_abb)) + geom_line() +
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(limits = c(-1.4, 2.2)) +
  geom_rug(y = NA, alpha = .1) +
  geom_dl(alpha = .5, method = list("first.points", cex = .7)) +
  labs(x = "Standardized average for 101 countries (k_c^')",
       y = "Expected prevalence of deprivation")
ggsave(paste0(print.images, "02_te_eff.pdf"), width = 6.5, height = 3)

# FIGURE 3 ----
# function to transform mean & sample size to a-b param of beta dist.
beta_par <- function (m, n) c(m * n, (1 - m) * n)
# mean .5 and sample size 100 should return a = 50, b = 50
beta_par(.5, 100)
# mean .3 and sample size .5 should return a = .15, b = .35
beta_par(.3, .5)
bpx <- seq(1e-2, 1 - 1e-2, 1e-4)

(data.frame(bpx = bpx, m5n100 = dbeta(bpx, beta_par(.5, 100)[1], beta_par(.5, 100)[2]),
            m5n5 = dbeta(bpx, beta_par(.5, 5)[1], beta_par(.5, 5)[2]),
            m5n.5 = dbeta(bpx, beta_par(.5, .5)[1], beta_par(.5, .5)[2])) %>%
    ggplot(aes(bpx, m5n100)) + geom_line() +
    geom_line(aes(y = m5n.5), linetype = 2) +
    geom_vline(xintercept = .5, alpha = .5, linetype = 2) +
    annotate("text", x = .55, y = 6, label = "Beta(.5, 100)", hjust = 0) +
    annotate("text", x = .025, y = 4, label = "Beta(.5, .5)", hjust = 0) +
    scale_x_continuous(labels = percent_format(), name = "mean") +
    scale_y_continuous(name = "density") +
    theme_classic()) +
  (data.frame(bpx = bpx, m1n100 = dbeta(bpx, beta_par(.3, 100)[1], beta_par(.3, 100)[2]),
              m1n3 = dbeta(bpx, beta_par(.3, 3)[1], beta_par(.3, 3)[2])) %>%
     ggplot(aes(bpx, m1n100)) + geom_line() +
     geom_line(aes(y = m1n3), linetype = 2) +
     geom_vline(xintercept = .3, alpha = .5, linetype = 2) +
     annotate("text", x = .35, y = 7.5, label = "Beta(.3, 100)", hjust = 0) +
     annotate("text", x = .5, y = 1.2, label = "Beta(.3, 3)", hjust = 0) +
     scale_x_continuous(labels = percent_format(), name = "mean") +
     theme_classic() +
     scale_y_continuous(name = "density"))
ggsave(paste0(print.images, "03_beta_dist_ex.pdf"), width = 6.5, height = 2.5)

# tau-equivalent model (beta dist.) ----
# Rerun tau-equiv model using within 0-1 outcome
(te.reg.01 <- glmmTMB(
  response.01 ~ 0 + item_t + (1 | id), X,
  dispformula = ~ 0 + item_t))
# Run model with beta_family()
(te.reg.beta <- glmmTMB(
  response.01 ~ 0 + item_t + (1 | id), X, beta_family(),
  dispformula = ~ 0 + item_t))
(glmmTMB(
  response.01 ~ 0 + item_t + (1 | id) + (1 | c(1:603)), X, beta_family(),
  dispformula = ~ 0 + item_t))

# Indicators means from both models
# Note inv. logit transformation for beta means:
round(cbind(
  Gaussian = fixef(te.reg.01)$cond,
  beta = plogis(fixef(te.reg.beta)$cond)), 3)

# Compare models heteroskedastic Gaussian & beta models ----
anova(te.reg.01, te.reg.beta)

# FIGURE 4 ----
# Extract random effects and standardize
X$ave.te.b <- (ranef(te.reg.beta)$cond$id[, 1])[X$id] / attr(VarCorr(te.reg.beta)$cond$id, "stddev")

data.frame(name = X$item_abb, ave = X$ave.te.b, fitted_unc = predict(te.reg.beta),
           fitted_con = fitted(te.reg.beta)) %>%
  pivot_longer(3:4, names_pattern = "(.*)_(.*)", names_to = c(".value", "cons")) %>%
  mutate(cons = ifelse(cons == "con", "Probabilities", "Unconstrained prediction")) %>%
  ggplot(aes(ave, fitted, group = name, label = name)) + geom_line() + facet_wrap(~ cons, scales = "free") +
  scale_x_continuous(limits = c(-1.6, 2.5)) +
  geom_dl(alpha = .5, method = list("last.points", cex = .7)) +
  # geom_dl(method = "last.points", alpha = .5) +
  labs(x = "Standardized average for 101 countries (k_c^')",
       y = "Expected deprivation") +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "04_te_eff_beta.pdf"), width = 6.5, height = 3)

# FIGURE 5 ----
(data.frame(x = seq(0, .5, 1e-4),
            density = dlnorm(seq(0, .5, 1e-4), log(.1), log(.25 / .1) / qnorm(.95))) %>%
   ggplot(aes(x, density)) + geom_line() + theme_classic() +
   scale_x_continuous(breaks = c(0, .1, .25, .4)) +
   geom_vline(xintercept = c(.1, .25), linetype = 2, alpha = .5) +
   theme(axis.title = element_blank()) +
   labs(tag = "A", subtitle = "Density")) +
  (data.frame(x = seq(0, .5, 1e-4),
              cum_dens = plnorm(seq(0, .5, 1e-4), log(.1), log(.25 / .1) / qnorm(.95))) %>%
     ggplot(aes(x, cum_dens)) + geom_line() + theme_classic() +
     scale_x_continuous(breaks = c(0, .1, .25, .4)) +
     scale_y_continuous(labels = percent_format(), breaks = c(.0, .25, .5, .75, .95)) +
     theme(axis.title.x = element_text(hjust = 1.5), axis.title.y = element_blank()) +
     geom_hline(yintercept = c(.5, .95), linetype = 2, alpha = .5) +
     geom_vline(xintercept = c(.1, .25), linetype = 2, alpha = .5) +
     labs(x = "A-priori expectations for lambda_i parameter",
          subtitle = "Cumulative distribution", tag = "B"))
ggsave(paste0(print.images, "05_log_normal_lambda.pdf"), width = 5, height = 2.5)

# CONGENERIC CFA ----
# Congeneric CFA in lavaan ----
# Creating error SD as square root of error variance
writeLines(cong.lav.form <- paste(
  paste("F =~ ", paste0(colnames(dat)[15:20], collapse = " + ")),
  "x1 ~~ a * x1\nx2 ~~ b * x2\nx3 ~~ c * x3\nx4 ~~ d * x4\nx5 ~~ e * x5\nx6 ~~ f * x6",
  "a1 := sqrt(a)\nb1 := sqrt(b)\nc1 := sqrt(c)\nd1 := sqrt(d)\ne1 := sqrt(e)\nf1 := sqrt(f)",
  sep = "\n"))
summary(cong.cfa <- cfa(
  cong.lav.form, dat, std.lv = TRUE, meanstructure = TRUE, missing = "ML"),
  fit.measures = TRUE)

# Congeneric regression in Stan ----
# Data list with references to equation 8
data <- list(
  alpha_scale = 1,  # scale of a parameter
  beta_scale = sqrt(.25) / qnorm(.975),  # scale of sigma_t
  sigma_scale = .5 / qnorm(.975) / qnorm(.975),  # scale of sigma_i
  lambda_median = log(.1),  # log-median of lambda_i
  lambda_scale = log(.25 / .1) / qnorm(.95),  # scale of lambda_i
  # item id number long-form, ids be contiguous from 1 to number of items
  item_id = X$item,
  Ni = max(X$item),  # number of items
  # respondent id number long-form, ids be contiguous from 1 to number of respondents
  resp_id = X$id,
  Np = max(X$id),  # number of respondents
  y = X$response.01,  # response with all values between 0 and 1 exclusive
  N = length(X$response),  # rows in data (long form)
  ret_yhat = 1,  # return predicted outcomes (1: yes, 0: no)
  ret_ll = 0  # return log-likelihood (1: yes, 0 = no)
)

# Compile Stan model (takes time)
cong.mod.ln <- cmdstan_model(file.path(paste0(stan.scripts, "cong_ln.stan")))
cong.fit.ln <- cong.mod.ln$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems
)
cong.fit.ln$cmdstan_diagnose()  # Run quick diagnostic checks
cong.fit.ln.rs <- read_stan_csv(cong.fit.ln$output_files())  # convert to Rstan

# Print major parameters
print(cong.fit.ln.rs, c("theta_p", "yhat", "lp__"), include = FALSE)
#            mean se_mean   sd 2.5%  25%  50%  75% 97.5% n_eff Rhat
# alpha      0.34       0 0.07 0.21 0.30 0.34 0.38  0.48  7160 1.00
# sigma_beta 0.15       0 0.06 0.08 0.11 0.14 0.17  0.30 16216 1.00
# beta[1]    0.51       0 0.03 0.44 0.49 0.51 0.53  0.58  1481 1.01
# beta[2]    0.42       0 0.03 0.36 0.40 0.42 0.44  0.47  1486 1.01
# beta[3]    0.27       0 0.02 0.23 0.25 0.27 0.28  0.30  1565 1.01
# beta[4]    0.31       0 0.03 0.25 0.29 0.31 0.33  0.37  1259 1.01
# beta[5]    0.39       0 0.03 0.33 0.37 0.39 0.41  0.44  1571 1.01
# beta[6]    0.17       0 0.02 0.14 0.16 0.17 0.19  0.21  1560 1.01
# lambda[1]  0.29       0 0.02 0.25 0.28 0.29 0.31  0.35  3913 1.00
# lambda[2]  0.25       0 0.02 0.21 0.23 0.25 0.26  0.29  3958 1.00
# lambda[3]  0.16       0 0.01 0.13 0.15 0.16 0.17  0.19  4034 1.00
# lambda[4]  0.29       0 0.02 0.25 0.27 0.29 0.30  0.33  2656 1.00
# lambda[5]  0.23       0 0.02 0.19 0.22 0.23 0.25  0.28  4229 1.00
# lambda[6]  0.15       0 0.01 0.13 0.14 0.15 0.16  0.18  3905 1.00
# sigma[1]   0.17       0 0.01 0.14 0.16 0.17 0.18  0.20 25203 1.00
# sigma[2]   0.14       0 0.01 0.12 0.13 0.14 0.15  0.17 26849 1.00
# sigma[3]   0.10       0 0.01 0.09 0.10 0.10 0.11  0.12 29277 1.00
# sigma[4]   0.09       0 0.01 0.07 0.08 0.09 0.10  0.11  7390 1.00
# sigma[5]   0.16       0 0.01 0.13 0.15 0.15 0.16  0.18 24396 1.00
# sigma[6]   0.09       0 0.01 0.08 0.09 0.09 0.10  0.11 29276 1.00

# Traceplots:
mcmc_trace(cong.fit.ln.rs, regex_pars = c("alpha", "beta", "lambda", "sigma"))
# Rank plots:
mcmc_rank_overlay(
  cong.fit.ln.rs, regex_pars = c("alpha", "beta", "lambda", "sigma"),
  ref_line = TRUE)
# Summarize draws
cong.pars.df <- as.data.frame(
  summary(cong.fit.ln.rs, c("beta", "lambda", "sigma"))$summary)
# beta[1-6] are means
# lambda[1-6] are loadings
# sigma[1-6] are residual standard deviations
cong.pars.df$variable <- rownames(cong.pars.df)

# FIGURE 6 ----
# Join lavaan and Bayes results
cong.pars.df.joined <- cbind(
  cong.pars.df[, c(11, 1, 6, 4, 8)],
  parameterEstimates(cong.cfa)[c(14:19, 1:6, 21:26), c(5, 9:10)])
cong.pars.df.joined
# Rename columns
names(cong.pars.df.joined) <- c(
  "variable", "mean_Bayesian-multilevel", "median", "lower_Bayesian-multilevel",
  "upper_Bayesian-multilevel", "mean_Lavaan", "lower_Lavaan", "upper_Lavaan")
cong.pars.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
cong.pars.df.joined %>%
  pivot_longer(c(2, 4:8), names_pattern = "(.*)_(.*)", names_to = c(".value", "method")) %>%
  mutate(type = gsub("\\[[0-9]\\]", "", variable), number = gsub("[a-z]+|\\[|\\]", "", variable)) %>%
  mutate(type = factor(type, levels = c("lambda", "beta", "sigma"),
                       labels = c("Loadings", "Prevalence", "Residual SD"))) %>%
  ggplot(aes(number, mean, fill = method, shape = method)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(.5), alpha = .65) +
  theme_classic() + theme(legend.position = "top", strip.background = element_blank(),
                          axis.title.y = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ type, scales = "free_y") +
  labs(
    x = "Indicator number", subtitle = "Estimates with 95% interval", shape = "Method", fill = "Method",
    caption = paste0(fig.capt, "\n",
                     "Quantile interval for Bayesian model, and normal-theory intervals for lavaan model"))
ggsave(paste0(print.images, "06_cong_params_lav_bayes.pdf"), width = 6.5, height = 4)

# FIGURE 7 ----
cong.pars.df.2 <- as.data.frame(summary(cong.fit.ln.rs, c("beta", "lambda", "theta_p"))$summary)
cong.pars.df.2 <- cong.pars.df.2[, c(1, 6)]
cong.pars.df.2$variable <- rownames(cong.pars.df.2)
cong.pars.df.2 <- as.data.table(cong.pars.df.2)
# Use equation to create predicted outcome in wide form
cong.pars.df.2.res <-
  t(sapply(cong.pars.df.2[variable %in% paste0("theta_p[", 1:101, "]"), mean], function (theta) {
    theta * cong.pars.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean] +
      cong.pars.df.2[variable %in% paste0("beta[", 1:6, "]"), mean]
  }))
colnames(cong.pars.df.2.res) <- levels(X$item_n)
# Add standardized averages and row id
cong.pars.df.2.res <- data.frame(
  cong.pars.df.2.res,
  theta = cong.pars.df.2[variable %in% paste0("theta_p[", 1:101, "]"), mean], id = 1:101)
# Transform data to long form
cong.pars.df.2.res <- cong.pars.df.2.res %>% pivot_longer(1:6) %>% as.data.frame()
cong.pars.df.2.res$item <- 1:6  # Add item IDs
# Add lambda parameter
cong.pars.df.2.res$lambda <- cong.pars.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean]
cong.pars.df.2.res %>%
  ggplot(aes(theta, value, group = name, label = item)) +
  geom_line(aes(alpha = lambda)) + geom_rug(y = NA, alpha = .1) +
  geom_dl(aes(alpha = lambda), method = "last.points",
          data = filter(cong.pars.df.2.res, item %in% c(1, 2, 3, 5, 6))) +
  geom_dl(aes(alpha = lambda), method = "first.points",
          data = filter(cong.pars.df.2.res, !item %in% c(1, 2, 3, 6))) +
  scale_alpha(range = c(.3, 1)) + guides(alpha = FALSE) +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 4, .2)) +
  theme(panel.grid.minor.y = element_blank()) +
  labs(y = "Expected prevalence of deprivation by indicator",
       x = "Standardized average for 101 countries (k_c^')",
       caption = fig.capt)
ggsave(paste0(print.images, "07_cong_eff.pdf"), width = 6.5, height = 4)

# FIGURE 8 ----
scenarios.df <- expand.grid(
  prev = seq(.05, .5, length.out = 4),  # Expected prevalences (probability)
  country = c(-2, 2),  # Country average on logit scale (two extremes)
  lambda = seq(0, 5, 1)  # Different values of loadings
)
# Create expected outcome on logit scale:
scenarios.df$lgt.exp <- with(scenarios.df, qlogis(prev) + country * lambda)
# Create expected outcome on probability scale:
scenarios.df$exp <- plogis(scenarios.df$lgt.exp)

scenarios.df %>%
  ggplot(aes(lambda, exp, group = country, linetype = factor(country))) + geom_point() + geom_line() +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ factor(percent(prev), percent(seq(.05, .5, length.out = 4))), nrow = 1) +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing = unit(.5, "cm"), axis.ticks = element_blank(),
        legend.position = "top", panel.grid.minor.x = element_blank()) +
  labs(x = "lambda (loading)", y = "Expected prevalence",
       linetype = "Value of standardized average",
       caption = "Each panel is a hypothetical indicator, and the labels (5%, 20%, ...) are the indicator means.",
       title = "Effect of loading on probabilities depends on indicator mean")
ggsave(paste0(print.images, "08_loading_logit.pdf"), width = 6.5, height = 4.5)

# Congeneric CFA in lavaan (bootstrapped) ----
summary(cong.cfa.b <- cfa(
  cong.lav.form, dat, std.lv = TRUE, meanstructure = TRUE,
  missing = "ML", se = "bootstrap",
  # Next line is to speed up bootstrap via parallelization
  bootstrap = 1e3, parallel = "multicore", ncpus = 6),
  fit.measures = TRUE, standardize = TRUE)

# Congeneric regression in Stan (beta dist.) ----
# Data list with references to equation 11
data.b <- list(
  alpha_scale = 5 / qnorm(.975),  # scale of a parameter
  beta_scale = 2 / qnorm(.975),  # scale of sigma_t
  lambda_median = log(.5),  # log-median of lambda_i
  lambda_scale = log(4 / .5) / qnorm(.975),  # scale of lambda_i
  # scaler element can be used to rescale the beta-probability to another scale
  # default is 1, but if scaler = 7, then response scale is assumed to be 0-7
  scaler = 1,
  # All other elements are same as Gaussian model above
  item_id = X$item, Ni = max(X$item), resp_id = X$id, Np = max(X$id),
  y = X$response.01, N = length(X$response), ret_yhat = 1, ret_ll = 0
)

# Compile Stan model
cong.mod.b <- cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_beta.stan")))
# Same configuration settings as Gaussian model
cong.fit.b <- cong.mod.b$sample(
  data = data.b, seed = 12345, iter_warmup = 1e3,
  iter_sampling = 3e3, chains = 8, parallel_chains = 8)

cong.fit.b$cmdstan_diagnose()
cong.fit.b.rs <- read_stan_csv(cong.fit.b$output_files())

print(cong.fit.b.rs, c("theta_p", "yhat", "lp__"), include = FALSE)
#             mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
# alpha      -0.91    0.01  0.40 -1.70 -1.16 -0.91 -0.67  -0.10  5234    1
# sigma_beta  0.88    0.00  0.29  0.48  0.67  0.82  1.03   1.61 18795    1
# beta[1]     0.03    0.01  0.19 -0.33 -0.10  0.03  0.15   0.39  1043    1
# beta[2]    -0.42    0.00  0.15 -0.72 -0.52 -0.42 -0.33  -0.13  1059    1
# beta[3]    -1.22    0.00  0.12 -1.45 -1.30 -1.22 -1.14  -0.99  1250    1
# beta[4]    -1.51    0.01  0.23 -1.96 -1.67 -1.51 -1.36  -1.06   885    1
# beta[5]    -0.53    0.00  0.14 -0.81 -0.62 -0.53 -0.43  -0.24  1203    1
# beta[6]    -1.94    0.00  0.15 -2.23 -2.04 -1.94 -1.84  -1.66  1189    1
# lambda[1]   1.65    0.00  0.16  1.36  1.54  1.65  1.76   1.99  3468    1
# lambda[2]   1.33    0.00  0.12  1.11  1.25  1.33  1.41   1.58  2979    1
# lambda[3]   0.95    0.00  0.10  0.77  0.88  0.94  1.01   1.14  3622    1
# lambda[4]   2.24    0.00  0.17  1.93  2.12  2.23  2.35   2.60  1964    1
# lambda[5]   1.17    0.00  0.12  0.95  1.08  1.16  1.24   1.41  3717    1
# lambda[6]   1.22    0.00  0.11  1.01  1.14  1.21  1.29   1.45  3043    1
# prec[1]     5.67    0.01  0.91  4.07  5.03  5.60  6.24   7.61 20916    1
# prec[2]    10.49    0.01  1.59  7.65  9.39 10.39 11.50  13.88 26728    1
# prec[3]    13.10    0.01  1.91  9.66 11.77 12.98 14.32  17.14 28339    1
# prec[4]    66.34    0.30 20.98 33.79 51.33 63.38 78.29 115.49  4997    1
# prec[5]     6.66    0.01  0.97  4.91  5.98  6.60  7.29   8.71 25933    1
# prec[6]    16.97    0.02  2.69 12.18 15.10 16.79 18.65  22.70 23136    1
# i_means[1]  0.51    0.00  0.05  0.42  0.48  0.51  0.54   0.60  1043    1
# i_means[2]  0.40    0.00  0.04  0.33  0.37  0.40  0.42   0.47  1053    1
# i_means[3]  0.23    0.00  0.02  0.19  0.21  0.23  0.24   0.27  1236    1
# i_means[4]  0.18    0.00  0.03  0.12  0.16  0.18  0.20   0.26   868    1
# i_means[5]  0.37    0.00  0.03  0.31  0.35  0.37  0.39   0.44  1195    1
# i_means[6]  0.13    0.00  0.02  0.10  0.12  0.13  0.14   0.16  1167    1
# sigma[1]    0.19    0.00  0.01  0.17  0.18  0.19  0.20   0.22 19314    1
# sigma[2]    0.14    0.00  0.01  0.13  0.14  0.14  0.15   0.17 13134    1
# sigma[3]    0.11    0.00  0.01  0.10  0.11  0.11  0.12   0.13  7149    1
# sigma[4]    0.05    0.00  0.01  0.03  0.04  0.05  0.05   0.07  2942    1
# sigma[5]    0.18    0.00  0.01  0.15  0.17  0.17  0.18   0.20 10820    1
# sigma[6]    0.08    0.00  0.01  0.06  0.07  0.08  0.08   0.10  3964    1

mcmc_trace(cong.fit.b.rs, regex_pars = c("i_means", "lambda", "sigma"))
mcmc_rank_overlay(
  cong.fit.b.rs, regex_pars = c("i_means", "lambda", "sigma"), ref_line = TRUE)

# Congeneric regression in Stan with OLRE (beta dist.) ----
cong.mod.b.olre <- cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_beta_olre.stan")))
# Included maxtreedepth settings for convergence
cong.fit.b.olre <- cong.mod.b.olre$sample(
  data = data.b, seed = 12345, iter_warmup = 1e3,
  iter_sampling = 3e3, chains = 8, parallel_chains = 8,
  max_treedepth = 15, adapt_delta = .99)

cong.fit.b.olre$cmdstan_diagnose()
# The E-BFMI, 0.122475, is below the nominal threshold of 0.3 which suggests that
# HMC may have trouble exploring the target distribution.
cong.fit.b.olre.rs <- read_stan_csv(cong.fit.b.olre$output_files())

print(cong.fit.b.olre.rs, c("theta_p", "yhat", "olre", "lp__"), include = FALSE, digits_summary = 4)
#                mean se_mean      sd    2.5%     25%      50%      75%    97.5% n_eff   Rhat
# alpha       -1.0525  0.0043  0.4399 -1.9136 -1.3310  -1.0628  -0.7790  -0.1581 10631 1.0003
# sigma_beta   0.9589  0.0023  0.3030  0.5371  0.7415   0.9019   1.1126   1.6996 17515 1.0000
# beta[1]      0.0612  0.0042  0.2526 -0.4253 -0.1062   0.0592   0.2278   0.5696  3641 1.0035
# beta[2]     -0.5639  0.0030  0.1701 -0.8983 -0.6782  -0.5644  -0.4507  -0.2268  3287 1.0036
# beta[3]     -1.3612  0.0021  0.1310 -1.6192 -1.4493  -1.3615  -1.2728  -1.1039  4027 1.0027
# beta[4]     -1.6094  0.0044  0.2328 -2.0658 -1.7660  -1.6109  -1.4529  -1.1523  2857 1.0043
# beta[5]     -0.6923  0.0030  0.1754 -1.0353 -0.8096  -0.6934  -0.5765  -0.3435  3393 1.0044
# beta[6]     -2.2885  0.0031  0.1719 -2.6265 -2.4036  -2.2867  -2.1736  -1.9521  3004 1.0036
# lambda[1]    2.2543  0.0026  0.2038  1.8901  2.1117   2.2434   2.3845   2.6851  6063 1.0009
# lambda[2]    1.5513  0.0018  0.1344  1.3054  1.4572   1.5458   1.6382   1.8320  5717 1.0011
# lambda[3]    1.0512  0.0014  0.1098  0.8492  0.9754   1.0460   1.1222   1.2792  6276 1.0009
# lambda[4]    2.3106  0.0026  0.1744  1.9955  2.1893   2.3012   2.4216   2.6806  4598 1.0014
# lambda[5]    1.4824  0.0018  0.1427  1.2226  1.3836   1.4755   1.5736   1.7801  6081 1.0008
# lambda[6]    1.5614  0.0022  0.1375  1.3126  1.4652   1.5554   1.6504   1.8499  3786 1.0013
# prec       111.4518  0.9508 28.1253 64.6501 91.4872 108.3350 128.3932 174.6634   875 1.0096
# nu          20.8983  0.0882 14.0566  3.0006 10.6050  17.7952  27.8835  56.0275 25372 1.0000
# tau          1.1261  0.0036  0.3920  0.5384  0.8715   1.0719   1.3076   2.0662 12064 1.0003
# olre_sd[1]   1.2528  0.0024  0.2367  0.8420  1.0942   1.2324   1.3894   1.7740 10038 1.0006
# olre_sd[2]   0.7690  0.0013  0.1464  0.5138  0.6711   0.7578   0.8527   1.0907 11929 1.0004
# olre_sd[3]   0.7836  0.0017  0.1506  0.5219  0.6823   0.7719   0.8698   1.1234  8180 1.0010
# olre_sd[4]   0.3235  0.0043  0.1553  0.0315  0.2167   0.3249   0.4273   0.6316  1322 1.0045
# olre_sd[5]   0.9736  0.0017  0.1830  0.6537  0.8503   0.9577   1.0794   1.3752 11096 1.0004
# olre_sd[6]   0.6243  0.0014  0.1269  0.4031  0.5382   0.6130   0.6985   0.9074  8069 1.0012
# i_means[1]   0.5151  0.0010  0.0621  0.3952  0.4735   0.5148   0.5567   0.6387  3636 1.0035
# i_means[2]   0.3635  0.0007  0.0391  0.2894  0.3367   0.3625   0.3892   0.4436  3272 1.0036
# i_means[3]   0.2049  0.0003  0.0213  0.1653  0.1901   0.2040   0.2188   0.2490  4009 1.0027
# i_means[4]   0.1692  0.0006  0.0327  0.1125  0.1460   0.1665   0.1895   0.2401  2831 1.0044
# i_means[5]   0.3346  0.0007  0.0388  0.2621  0.3080   0.3333   0.3597   0.4150  3403 1.0043
# i_means[6]   0.0931  0.0003  0.0145  0.0675  0.0829   0.0922   0.1021   0.1243  2959 1.0036

mcmc_trace(cong.fit.b.olre.rs, regex_pars = c("i_means", "lambda"))
mcmc_rank_overlay(
  cong.fit.b.olre.rs, regex_pars = c("i_means", "lambda"), ref_line = TRUE)

# FIGURE 9 ----
# Join lavaan and Bayes results
cong.pars.b.df <- as.data.frame(
  summary(cong.fit.b.olre.rs, c("i_means", "lambda", "olre_sd"))$summary)
# i_means[1-6] are means on probability scale
# lambda[1-6] are loadings
# olre_sd[1-6] are residual standard deviations
cong.pars.b.df$variable <- rownames(cong.pars.b.df)

cong.pars.b.df.joined <- cbind(
  cong.pars.b.df[, c(11, 1, 6, 4, 8)],
  parameterEstimates(cong.cfa.b)[c(14:19, 1:6, 21:26), c(5, 9:10)])
# Rename columns
names(cong.pars.b.df.joined) <-
  c("variable", "mean_Bayesian-multilevel-beta (OLRE)", "median", "lower_Bayesian-multilevel-beta (OLRE)",
    "upper_Bayesian-multilevel-beta (OLRE)", "mean_Lavaan", "lower_Lavaan", "upper_Lavaan")
cong.pars.b.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
cong.pars.b.df.joined %>%
  pivot_longer(c(2, 4:8), names_pattern = "(.*)_(.*)", names_to = c(".value", "method")) %>%
  mutate(type = gsub("\\[[0-9]\\]", "", variable), number = gsub("_|[a-z]+|\\[|\\]", "", variable)) %>%
  mutate(type = factor(type, levels = c("lambda", "i_means", "olre_sd"),
                       labels = c("Loadings", "Prevalence", "Residual SD"))) %>%
  mutate(type = case_when(
    type == "Loadings" & method == "Lavaan" ~ "Loadings (lavaan)",
    type == "Loadings" & method != "Lavaan" ~ "Loadings (Bayes)",
    type == "Residual SD" & method == "Lavaan" ~ "Residual SD (lavaan)",
    type == "Residual SD" & method != "Lavaan" ~ "Scale (Bayes)",
    TRUE ~ as.character(type)
  )) %>%
  mutate(type = factor(type, c("Loadings (Bayes)", "Loadings (lavaan)", "Prevalence",
                               "Scale (Bayes)", "Residual SD (lavaan)"))) %>%
  ggplot(aes(number, mean, fill = method, shape = method)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(.5), alpha = .65) +
  theme_classic() + theme(legend.position = "top", strip.background = element_blank(),
                          axis.title.y = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ type, scales = "free_y") +
  labs(
    x = "Indicator number", subtitle = "Estimates with 95% interval", shape = "Method", fill = "Method",
    caption = paste0(fig.capt, "\n",
                     "Quantile interval for Bayesian model, and bootstrap intervals for lavaan model"))
ggsave(paste0(print.images, "09_cong_params_lav_bayes_beta.pdf"), width = 6.5, height = 5.5)

# FIGURE 10 ----
# Extract predicted outcomes: 100 samples randomly drawn from 24000 samples
G.cong.beta.olre <- as.data.frame(cong.fit.b.olre.rs, "yhat")
# Should be 100 samples by 603 responses
dim(G.cong.beta.olre <- G.cong.beta.olre[sample(1:nrow(G.cong.beta.olre), 100), ])
G.cong.beta.olre <- data.frame(item = X$item_n, response = X$response.01, t(G.cong.beta.olre)) %>%
  pivot_longer(contains("X"), names_pattern = "X(.*)")

G.cong.beta.olre %>%
  ggplot(aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(labels = percent_format(), name = "Distribution of prevalence") +
  geom_density(aes(value, group = name), fill = NA, col = "darkgrey") +
  geom_density(col = 1, fill = NA, data = filter(G.cong.beta.olre, name == min(as.integer(name)))) +
  theme_classic() +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "10_cong_eff_beta_ppc.pdf"), width = 6.5, height = 4)

# FIGURE 11 ----
# Extract item means
i.means.df <- as.data.frame(cong.fit.b.olre.rs, "i_means")
# Rank each row
tmp.res <- apply(i.means.df, 1, rank)
(tmp.res <- as.data.frame(t(apply(tmp.res, 1, function (row) {
  sapply(1:max(tmp.res), function (i) sum(row == i))
}))))
tmp.res$indicator <- factor(levels(X$item_n), levels(X$item_n))
tmp.res

plt.a <- tmp.res %>% pivot_longer(1:6) %>%
  group_by(indicator) %>%
  mutate(value = ifelse(is.na(value), 0, value),
         value = value / sum(value),
         value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(indicator, name, label = percent(value, .1))) + geom_tile(fill = NA, col = "lightgrey") +
  scale_y_discrete(labels = 1:6) +
  scale_alpha(range = c(.25, 1)) + theme_classic() + geom_text(aes(alpha = value)) +
  guides(alpha = FALSE) + theme(axis.title.x = element_blank()) +
  labs(subtitle = "Probalistic ranking of indicators (1 = lowest, 6 = highest prevalence of deprivation)",
       tag = "A", y = "Lowest to highest")
plt.a

# Extract loadings following same set-up as above:
loads.df <- as.data.frame(cong.fit.b.olre.rs, "lambda")
tmp.res <- apply(loads.df, 1, rank)
(tmp.res <- as.data.frame(t(apply(tmp.res, 1, function (row) {
  sapply(1:max(tmp.res), function (i) sum(row == i))
}))))
tmp.res$indicator <- factor(levels(X$item_n), levels(X$item_n))
tmp.res

plt.b <- tmp.res %>% pivot_longer(1:6) %>%
  group_by(indicator) %>%
  mutate(value = ifelse(is.na(value), 0, value),
         value = value / sum(value),
         value = ifelse(value == 0, NA, value)) %>%
  ggplot(aes(indicator, name, label = percent(value, .1))) + geom_tile(fill = NA, col = "lightgrey") +
  scale_alpha(range = c(.25, 1)) + theme_classic() + geom_text(aes(alpha = value)) +
  scale_y_discrete(labels = 1:6) +
  guides(alpha = FALSE) +
  labs(subtitle = "Probalistic ranking of indicators (1 = weakest, 6 = strongest relation to living standards)",
       tag = "B", y = "Weakest to strongest", x = "Indicator")
plt.a / plt.b
ggsave(paste0(print.images, "11_rank_indicators.pdf"), width = 6.5, height = 4)

# FIGURE 12 ----
# follows same approach as in FIGURE 7
cong.pars.b.df.2 <- as.data.frame(summary(cong.fit.b.olre.rs, c("beta", "lambda", "theta_p"))$summary)
cong.pars.b.df.2 <- cong.pars.b.df.2[, c(1, 6)]
cong.pars.b.df.2$variable <- rownames(cong.pars.b.df.2)
cong.pars.b.df.2 <- as.data.table(cong.pars.b.df.2)
cong.pars.b.df.2.res <-
  t(sapply(cong.pars.b.df.2[variable %in% paste0("theta_p[", 1:101, "]"), mean], function (theta) {
    theta * cong.pars.b.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean] +
      cong.pars.b.df.2[variable %in% paste0("beta[", 1:6, "]"), mean]
  }))
cong.pars.b.df.2.res <- data.frame(
  cong.pars.b.df.2.res,
  theta = cong.pars.b.df.2[variable %in% paste0("theta_p[", 1:101, "]"), mean], id = 1:101)
names(cong.pars.b.df.2.res)[1:6] <- levels(X$item_n)
cong.pars.b.df.2.res <- as.data.frame(cong.pars.b.df.2.res) %>% pivot_longer(1:6) %>% as.data.frame()
cong.pars.b.df.2.res$item <- 1:6
cong.pars.b.df.2.res$lambda <- cong.pars.b.df.2[variable %in% paste0("lambda[", 1:6, "]"), mean]
cong.pars.b.df.2.res$prob <- plogis(cong.pars.b.df.2.res$value)

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
cong.pars.b.df.2.res %>%
  ggplot(aes(theta * -1, prob, group = name, label = item)) +
  geom_line(aes(alpha = lambda, size = lambda)) + geom_rug(y = NA, alpha = .1) +
  geom_dl(aes(alpha = lambda), method = "first.points",
          data = filter(cong.pars.b.df.2.res, item %in% c(1, 2, 3, 5, 6))) +
  geom_dl(aes(alpha = lambda), method = "last.points", data = filter(cong.pars.b.df.2.res, item == 4)) +
  scale_alpha(range = c(.2, 1)) + scale_size(range = c(.25, 1)) +
  guides(alpha = FALSE, size = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme(panel.grid.minor.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Expected prevalence of deprivation by indicator",
       x = "Living standards (standardized averages multiplied by -1)",
       caption = fig.capt)
ggsave(paste0(print.images, "12_cong_eff_beta.pdf"), width = 6.5, height = 4)

# FIGURE 13 ----
# Extract country averages: should 24000 by 101 countries
dim(theta.c.b <- as.data.frame(cong.fit.b.olre.rs, "theta_p"))
# Calculate % each country has bottom 30 scores
ret.c.b <- rowMeans(apply(theta.c.b, 1, function (xs) {
  top.30 <- sort(xs, decreasing = TRUE)[1:30]
  sapply(xs, function (x) { x %in% top.30 })
}))
# Join bottom-30 % & country means
ret.c.b <- data.frame(ave = colMeans(theta.c.b), prob.30 = ret.c.b)
# Add bottom 30 indicator
ret.c.b$low.30 <- ret.c.b$ave %in% sort(ret.c.b$ave, decreasing = TRUE)[1:30]
# Add country ID
ret.c.b$ID <- 1:nrow(ret.c.b)

# Number of countries with >99% of bottom 30 score
sum(ret.c.b$prob.30 > .99)
# [1] 18
sum(ret.c.b$prob.30 > .95)  # 95%
# [1] 21
sum(ret.c.b$prob.30 > .90)  # 90%
# [1] 24
sum(ret.c.b$prob.30 > .50)  # 50%
# [1] 31

ret.c.b %>%
  ggplot(aes(-1 * ave, prob.30, shape = low.30, label = ID)) + geom_point(size = 3) +
  geom_hline(yintercept = c(.5, .8, .9, .95), size = .25, linetype = 2, alpha = .25) +
  geom_text_repel(data = filter(ret.c.b, prob.30 > .5 & prob.30 <= .95), nudge_x = .5, segment.alpha = .25) +
  geom_text_repel(aes(label = percent(prob.30, 1)),
                  data = filter(ret.c.b, prob.30 > .5 & prob.30 <= .8), nudge_x = -.2, segment.alpha = .25) +
  annotate("text", x = 1, y = .95, label = "95%") +
  scale_y_continuous(labels = percent_format(), breaks = c(0, .5, .8, .9, 1)) + theme_classic() +
  scale_shape_manual(values = c(4, 1)) +
  guides(shape = guide_legend(reverse = TRUE), linetype = FALSE) +
  theme(legend.position = "top", axis.title.y = element_blank()) +
  labs(x = "Living standards (standardized averages * -1)",
       subtitle = "Probability of having a bottom 30 living standard",
       shape = "Average score in bottom 30?")
ggsave(paste0(print.images, "13_bottom_30_cong_beta.pdf"), width = 6.5, height = 4)

# Comparing country 81 and 82 posterior averages
theta.9581 <- as.data.frame(cong.fit.b.olre.rs, c("theta_p[95]", "theta_p[81]"))
mean(theta.9581[, 1] > theta.9581[, 2])
# [1] 0.480875
theta.8382 <- as.data.frame(cong.fit.b.olre.rs, c("theta_p[83]", "theta_p[82]"))
mean(theta.8382[, 1] > theta.8382[, 2])  # 81 only marginally higher probability than 82
# [1] 0.5343333

# APPENDIX CONTENT ----

# Rankplots for congeneric model, Gaussian
mcmc_rank_overlay(cong.fit.ln.rs, regex_pars = c("alpha", "beta", "lambda", "sigma"),
                  ref_line = TRUE)
ggsave(paste0(print.images, "app_01_rank_gauss_cong.pdf"), height = 6, width = 8)
# Rankplots for congeneric model, Beta
mcmc_rank_overlay(cong.fit.b.olre.rs, regex_pars = c("i_means", "lambda", "olre_sd"),
                  ref_line = TRUE)
ggsave(paste0(print.images, "app_02_rank_beta_cong.pdf"), height = 6, width = 8)

# Posterior predictive for congeneric CFA as regression (Gaussian)
G.cong.ln <- as.data.frame(cong.fit.ln.rs, paste0("yhat"))
# Should be 100 samples by 603 responses
dim(G.cong.ln <- G.cong.ln[sample(1:nrow(G.cong.ln), 100), ])
G.cong.ln <- data.frame(item = X$item_n, response = X$response.01, t(G.cong.ln)) %>%
  pivot_longer(contains("X"), names_pattern = "X(.*)")

G.cong.ln %>%
  ggplot(aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(labels = percent_format(), name = "Distribution of prevalence") +
  geom_density(aes(value, group = name), fill = NA, col = "darkgrey") +
  geom_density(col = 1, fill = NA, data = filter(G.cong.ln, name == min(as.integer(name)))) +
  theme_classic() +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "app_03_ppc_gauss_cong.pdf"), height = 4, width = 6)

# Posterior predictive for congeneric beta CFA (NO OLRE) as regression (Equation 9)
# Extract predicted outcomes: 100 samples randomly drawn from 24000 samples
G.cong.beta <- as.data.frame(cong.fit.b.rs, "yhat")
# Should be 100 samples by 603 responses
dim(G.cong.beta <- G.cong.beta[sample(1:nrow(G.cong.beta), 100), ])
G.cong.beta <- data.frame(item = X$item_n, response = X$response.01, t(G.cong.beta)) %>%
  pivot_longer(contains("X"), names_pattern = "X(.*)")

G.cong.beta %>%
  ggplot(aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(labels = percent_format(), name = "Distribution of prevalence") +
  geom_density(aes(value, group = name), fill = NA, col = "darkgrey") +
  geom_density(col = 1, fill = NA, data = filter(G.cong.beta, name == min(as.integer(name)))) +
  theme_classic() +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "app_04_ppc_beta_cong.pdf"), width = 6.5, height = 4)

# More efficient estimation (of beta congeneric model) via "marginal likelihood" ----

# Compile Stan model
cong.mod.b.m <- cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_beta_marg.stan")))
# Same configuration settings as Gaussian model BUT 1K post-warmup samples / chain
cong.fit.b.m <- cong.mod.b.m$sample(
  data = data.b, seed = 12345, iter_warmup = 1e3,
  iter_sampling = 1e3, chains = 8, parallel_chains = 8)

cong.fit.b.m$cmdstan_diagnose()
cong.fit.b.m.rs <- read_stan_csv(cong.fit.b.m$output_files())

print(cong.fit.b.m.rs, c("latent_mv", "lp__"), include = FALSE)
#             mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
# alpha      -0.92    0.00  0.40 -1.70 -1.16 -0.92 -0.67  -0.09  7130    1
# sigma_beta  0.88    0.00  0.29  0.48  0.68  0.83  1.03   1.61  5362    1
# beta[1]     0.02    0.00  0.19 -0.34 -0.11  0.02  0.15   0.40  9958    1
# beta[2]    -0.43    0.00  0.15 -0.73 -0.53 -0.43 -0.33  -0.14 11038    1
# beta[3]    -1.22    0.00  0.12 -1.47 -1.30 -1.22 -1.14  -1.00  7947    1
# beta[4]    -1.52    0.00  0.23 -1.97 -1.68 -1.52 -1.37  -1.08 11210    1
# beta[5]    -0.53    0.00  0.14 -0.81 -0.63 -0.53 -0.44  -0.25  7745    1
# beta[6]    -1.96    0.00  0.15 -2.25 -2.05 -1.96 -1.86  -1.68  8103    1
# lambda[1]   1.66    0.00  0.17  1.35  1.55  1.65  1.77   2.01  5761    1
# lambda[2]   1.34    0.00  0.12  1.11  1.25  1.33  1.42   1.61  7767    1
# lambda[3]   0.95    0.00  0.10  0.77  0.88  0.94  1.01   1.15  6173    1
# lambda[4]   2.24    0.00  0.18  1.93  2.12  2.24  2.36   2.62  9185    1
# lambda[5]   1.17    0.00  0.12  0.95  1.09  1.17  1.25   1.43  6435    1
# lambda[6]   1.23    0.00  0.12  1.02  1.15  1.22  1.30   1.47  7277    1
# prec[1]     5.73    0.01  0.93  4.08  5.08  5.68  6.31   7.71  5582    1
# prec[2]    10.81    0.02  1.66  7.93  9.63 10.71 11.85  14.39  8681    1
# prec[3]    13.41    0.02  2.01  9.77 12.00 13.29 14.68  17.62  7202    1
# prec[4]    68.92    0.45 21.39 36.17 53.29 65.83 81.40 118.41  2216    1
# prec[5]     6.81    0.01  1.00  5.02  6.12  6.75  7.45   8.94  6729    1
# prec[6]    17.78    0.04  2.91 12.62 15.71 17.57 19.67  23.84  6437    1
# i_means[1]  0.51    0.00  0.05  0.42  0.47  0.51  0.54   0.60  9965    1
# i_means[2]  0.39    0.00  0.04  0.32  0.37  0.39  0.42   0.46 11048    1
# i_means[3]  0.23    0.00  0.02  0.19  0.21  0.23  0.24   0.27  7960    1
# i_means[4]  0.18    0.00  0.03  0.12  0.16  0.18  0.20   0.25 11034    1
# i_means[5]  0.37    0.00  0.03  0.31  0.35  0.37  0.39   0.44  7724    1
# i_means[6]  0.12    0.00  0.02  0.10  0.11  0.12  0.13   0.16  8008    1
# sigma[1]    0.19    0.00  0.01  0.17  0.18  0.19  0.20   0.22  5297    1
# sigma[2]    0.14    0.00  0.01  0.12  0.14  0.14  0.15   0.16  8819    1
# sigma[3]    0.11    0.00  0.01  0.09  0.10  0.11  0.12   0.13  6611    1
# sigma[4]    0.05    0.00  0.01  0.03  0.04  0.05  0.05   0.07  2583    1
# sigma[5]    0.17    0.00  0.01  0.15  0.17  0.17  0.18   0.20  6474    1
# sigma[6]    0.08    0.00  0.01  0.06  0.07  0.08  0.08   0.09  6042    1

mcmc_rank_overlay(cong.fit.b.m.rs, regex_pars = c("i_means", "lambda", "sigma"),
                  ref_line = TRUE)

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
