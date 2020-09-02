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
ggsave(paste0(print.images, "01_deprivation.pdf"), width = 6, height = 4)

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
  geom_dl(method = "first.points", alpha = .5) +
  labs(x = "Standardized average for 101 countries (k_c^')",
       y = "Expected prevalence of deprivation by indicator")
ggsave(paste0(print.images, "02_te_eff.pdf"), width = 6.5, height = 4)

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
ggsave(paste0(print.images, "03_beta_dist_ex.pdf"), width = 6.5, height = 4)

# tau-equivalent model (beta dist.) ----
# Rerun tau-equiv model using within 0-1 outcome
(te.reg.01 <- glmmTMB(
  response.01 ~ 0 + item_t + (1 | id), X,
  dispformula = ~ 0 + item_t))
# Run model with beta_family()
(te.reg.beta <- glmmTMB(
  response.01 ~ 0 + item_t + (1 | id), X, beta_family(),
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
  geom_dl(method = "last.points", alpha = .5) +
  labs(x = "Standardized average for 101 countries (k_c^')",
       y = "Expected deprivation by indicator") +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "04_te_eff_beta.pdf"), width = 6.5, height = 4)

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
# saveRDS(cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_sz.stan"))),
#         paste0(stan.scripts, "cong_ln_sz.rds"))
# (cong.mod.ln.sz <- readRDS(paste0(stan.scripts, "cong_ln_sz.rds")))
saveRDS(cmdstan_model(file.path(paste0(stan.scripts, "cong_ln.stan"))),
        paste0(stan.scripts, "cong_ln.rds"))
(cong.mod.ln <- readRDS(paste0(stan.scripts, "cong_ln.rds")))
cong.fit.ln <- cong.mod.ln$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  num_warmup = 1e3,  # 1000 warmup samples per chain
  num_samples = 3e3,  # 3000 samples for inference per chain
  num_chains = 8,  # 8 chains
  num_cores = 8  # 8 cores on multicore systems
)
cong.fit.ln$cmdstan_diagnose()  # Run quick diagnostic checks
cong.fit.ln.rstan <- read_stan_csv(cong.fit.ln$output_files())  # convert to Rstan

# Print major parameters
print(cong.fit.ln.rstan, c("theta_p", "yhat"), include = FALSE)
#              mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# alpha        0.34     0.0  0.07   0.21   0.30   0.34   0.38   0.48  8192 1.00
# sigma_beta   0.15     0.0  0.06   0.08   0.11   0.14   0.18   0.31 15903 1.00
# beta[1]      0.51     0.0  0.03   0.44   0.49   0.51   0.53   0.58  1592 1.01
# beta[2]      0.42     0.0  0.03   0.36   0.40   0.42   0.43   0.47  1518 1.01
# beta[3]      0.27     0.0  0.02   0.23   0.25   0.27   0.28   0.30  1642 1.01
# beta[4]      0.31     0.0  0.03   0.25   0.29   0.31   0.33   0.37  1329 1.01
# beta[5]      0.39     0.0  0.03   0.33   0.37   0.39   0.40   0.44  1743 1.01
# beta[6]      0.17     0.0  0.02   0.14   0.16   0.17   0.19   0.21  1623 1.01
# lambda[1]    0.30     0.0  0.02   0.25   0.28   0.29   0.31   0.35  4194 1.00
# lambda[2]    0.25     0.0  0.02   0.21   0.23   0.25   0.26   0.29  4002 1.00
# lambda[3]    0.16     0.0  0.01   0.13   0.15   0.16   0.17   0.19  4214 1.00
# lambda[4]    0.29     0.0  0.02   0.25   0.28   0.29   0.30   0.33  2848 1.00
# lambda[5]    0.23     0.0  0.02   0.19   0.22   0.23   0.25   0.28  4515 1.00
# lambda[6]    0.15     0.0  0.01   0.13   0.14   0.15   0.16   0.18  4210 1.00
# sigma[1]     0.17     0.0  0.01   0.14   0.16   0.17   0.18   0.20 23340 1.00
# sigma[2]     0.14     0.0  0.01   0.12   0.13   0.14   0.15   0.17 24093 1.00
# sigma[3]     0.10     0.0  0.01   0.09   0.10   0.10   0.11   0.12 26840 1.00
# sigma[4]     0.09     0.0  0.01   0.07   0.08   0.09   0.10   0.11  6453 1.00
# sigma[5]     0.16     0.0  0.01   0.13   0.15   0.15   0.16   0.18 21259 1.00
# sigma[6]     0.09     0.0  0.01   0.08   0.09   0.09   0.10   0.11 26331 1.00

# Traceplots:
mcmc_trace(cong.fit.ln.rstan, regex_pars = c("alpha", "beta", "lambda", "sigma"))
# Rank plots:
mcmc_rank_overlay(
  cong.fit.ln.rstan, regex_pars = c("alpha", "beta", "lambda", "sigma"),
  ref_line = TRUE)
# Summarize draws
cong.pars.df <- as.data.frame(
  summary(cong.fit.ln.rstan, c("beta", "lambda", "sigma"))$summary)
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
cong.pars.df.2 <- as.data.frame(summary(cong.fit.ln.rstan, c("beta", "lambda", "theta_p"))$summary)
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
# Data list with references to equation 10
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
saveRDS(cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_beta_sample_size.stan"))),
        paste0(stan.scripts, "cong_ln_beta_sample_size.rds"))
(cong.mod.b <- readRDS(paste0(stan.scripts, "cong_ln_beta_sample_size.rds")))
# Same configuration settings as Gaussian model
cong.fit.b <- cong.mod.b$sample(
  data = data.b, seed = 12345, num_warmup = 1e3,
  num_samples = 3e3, num_chains = 8, num_cores = 8)

cong.fit.b$cmdstan_diagnose()
cong.fit.b.rstan <- read_stan_csv(cong.fit.b$output_files())

print(cong.fit.b.rstan, c("theta_p", "yhat"), include = FALSE)
#              mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# alpha       -0.92    0.01  0.40  -1.72  -1.17  -0.92  -0.68  -0.12  4421 1.00
# sigma_beta   0.87    0.00  0.29   0.48   0.67   0.82   1.02   1.61 17545 1.00
# beta[1]      0.02    0.01  0.19  -0.35  -0.11   0.02   0.14   0.38   921 1.02
# beta[2]     -0.44    0.01  0.15  -0.73  -0.54  -0.43  -0.33  -0.15   851 1.02
# beta[3]     -1.23    0.00  0.12  -1.46  -1.30  -1.22  -1.15  -1.00  1047 1.01
# beta[4]     -1.53    0.01  0.23  -1.99  -1.69  -1.53  -1.37  -1.09   738 1.02
# beta[5]     -0.54    0.00  0.14  -0.82  -0.63  -0.54  -0.44  -0.26  1018 1.01
# beta[6]     -1.95    0.00  0.15  -2.24  -2.05  -1.95  -1.85  -1.67   991 1.01
# lambda[1]    1.65    0.00  0.16   1.35   1.54   1.65   1.76   1.99  3411 1.00
# lambda[2]    1.33    0.00  0.12   1.11   1.25   1.33   1.41   1.58  3010 1.00
# lambda[3]    0.94    0.00  0.10   0.77   0.88   0.94   1.01   1.14  3483 1.00
# lambda[4]    2.24    0.00  0.17   1.92   2.12   2.23   2.35   2.60  2046 1.00
# lambda[5]    1.16    0.00  0.12   0.94   1.08   1.16   1.24   1.41  3733 1.00
# lambda[6]    1.22    0.00  0.11   1.01   1.14   1.21   1.29   1.45  2981 1.00
# prec[1]      5.67    0.01  0.91   4.03   5.04   5.62   6.25   7.62 18392 1.00
# prec[2]     10.50    0.01  1.55   7.72   9.41  10.40  11.49  13.80 23372 1.00
# prec[3]     13.08    0.01  1.89   9.64  11.78  12.98  14.28  17.07 21987 1.00
# prec[4]     66.74    0.30 21.16  33.67  51.47  63.80  79.11 115.77  4898 1.00
# prec[5]      6.66    0.01  0.96   4.91   5.99   6.59   7.26   8.71 22086 1.00
# prec[6]     16.93    0.02  2.72  12.15  15.01  16.76  18.62  22.79 20913 1.00
# i_means[1]   0.50    0.00  0.05   0.41   0.47   0.50   0.54   0.59   921 1.02
# i_means[2]   0.39    0.00  0.04   0.32   0.37   0.39   0.42   0.46   845 1.02
# i_means[3]   0.23    0.00  0.02   0.19   0.21   0.23   0.24   0.27  1035 1.01
# i_means[4]   0.18    0.00  0.03   0.12   0.16   0.18   0.20   0.25   722 1.02
# i_means[5]   0.37    0.00  0.03   0.31   0.35   0.37   0.39   0.44  1008 1.01
# i_means[6]   0.13    0.00  0.02   0.10   0.11   0.12   0.14   0.16   965 1.01
# sigma[1]     0.19    0.00  0.01   0.17   0.18   0.19   0.20   0.22 16957 1.00
# sigma[2]     0.14    0.00  0.01   0.13   0.14   0.14   0.15   0.17  9882 1.00
# sigma[3]     0.11    0.00  0.01   0.10   0.11   0.11   0.12   0.13  5165 1.00
# sigma[4]     0.05    0.00  0.01   0.03   0.04   0.05   0.05   0.07  2651 1.01
# sigma[5]     0.18    0.00  0.01   0.15   0.17   0.17   0.18   0.20  9887 1.00
# sigma[6]     0.08    0.00  0.01   0.06   0.07   0.08   0.08   0.10  3184 1.00

mcmc_trace(cong.fit.b.rstan, regex_pars = c("i_means", "lambda", "sigma"))
mcmc_rank_overlay(
  cong.fit.b.rstan, regex_pars = c("i_means", "lambda", "sigma"), ref_line = TRUE)

cong.pars.b.df <- as.data.frame(
  summary(cong.fit.b.rstan, c("i_means", "lambda", "sigma"))$summary)
# i_means[1-6] are means on probability scale
# lambda[1-6] are loadings
# sigma[1-6] are residual standard deviations (probability scale)
cong.pars.b.df$variable <- rownames(cong.pars.b.df)

# FIGURE 9 ----
# Join lavaan and Bayes results
cong.pars.b.df.joined <- cbind(
  cong.pars.b.df[, c(11, 1, 6, 4, 8)],
  parameterEstimates(cong.cfa.b)[c(14:19, 1:6, 21:26), c(5, 9:10)])
# Rename columns
names(cong.pars.b.df.joined) <-
  c("variable", "mean_Bayesian-multilevel-beta", "median", "lower_Bayesian-multilevel-beta",
    "upper_Bayesian-multilevel-beta", "mean_Lavaan", "lower_Lavaan", "upper_Lavaan")
cong.pars.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
cong.pars.b.df.joined %>%
  pivot_longer(c(2, 4:8), names_pattern = "(.*)_(.*)", names_to = c(".value", "method")) %>%
  mutate(type = gsub("\\[[0-9]\\]", "", variable), number = gsub("_|[a-z]+|\\[|\\]", "", variable)) %>%
  mutate(type = factor(type, levels = c("lambda", "i_means", "sigma"),
                       labels = c("Loadings", "Prevalence", "Residual SD"))) %>%
  mutate(type = case_when(
    type == "Loadings" & method == "Lavaan" ~ "Loadings (lavaan)",
    type == "Loadings" & method != "Lavaan" ~ "Loadings (Bayes)",
    TRUE ~ as.character(type)
  )) %>%
  mutate(type = factor(type, c("Loadings (Bayes)", "Loadings (lavaan)", "Prevalence", "Residual SD"))) %>%
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
G.cong.beta <- as.data.frame(cong.fit.b.rstan, "yhat")
# Should be 100 samples by 603 responses
dim(G.cong.beta <- G.cong.beta[sample(1:nrow(G.cong.beta), 100), ])

data.frame(item = X$item_n, response = X$response.01, t(G.cong.beta)) %>%
  pivot_longer(contains("X"), names_pattern = "X(.*)") %>%
  ggplot(aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(labels = percent_format(), name = "Distribution of prevalence") +
  geom_density(aes(value, group = name), fill = NA, col = "darkgrey") +
  geom_density(col = 1, fill = NA) +
  theme_classic() +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "10_cong_eff_beta_ppc.pdf"), width = 6.5, height = 4)

# FIGURE 11 ----
# Extract item means
i.means.df <- as.data.frame(cong.fit.b.rstan, "i_means")
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
loads.df <- as.data.frame(cong.fit.b.rstan, "lambda")
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
cong.pars.b.df.2 <- as.data.frame(summary(cong.fit.b.rstan, c("beta", "lambda", "theta_p"))$summary)
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
  geom_line(aes(alpha = lambda)) + geom_rug(y = NA, alpha = .1) +
  geom_dl(aes(alpha = lambda), method = "first.points",
          data = filter(cong.pars.b.df.2.res, item %in% c(1, 2, 3, 5, 6))) +
  geom_dl(aes(alpha = lambda), method = "last.points", data = filter(cong.pars.b.df.2.res, item == 4)) +
  scale_alpha(range = c(.2, 1)) + guides(alpha = FALSE) +
  scale_y_continuous(labels = percent_format()) +
  theme(panel.grid.minor.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Expected prevalence of deprivation by indicator",
       x = "Living standards (standardized averages multiplied by -1)",
       caption = fig.capt)
ggsave(paste0(print.images, "12_cong_eff_beta.pdf"), width = 6.5, height = 4)

# FIGURE 13 ----
# Extract country averages: should 24000 by 101 countries
dim(theta.c.b <- as.data.frame(cong.fit.b.rstan, "theta_p"))
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
# [1] 20
sum(ret.c.b$prob.30 > .95)  # 95%
# [1] 24
sum(ret.c.b$prob.30 > .90)  # 90%
# [1] 26
sum(ret.c.b$prob.30 > .50)  # 50%
# [1] 31

ret.c.b %>%
  ggplot(aes(-1 * ave, prob.30, shape = low.30, label = ID)) + geom_point(size = 3) +
  geom_hline(yintercept = c(.5, .8, .9, .95), size = .25, linetype = 2, alpha = .25) +
  geom_text_repel(data = filter(ret.c.b, prob.30 > .5 & prob.30 <= .95), nudge_x = .5, segment.alpha = .25) +
  geom_text_repel(aes(label = percent(prob.30, 1)),
                  data = filter(ret.c.b, prob.30 > .5 & prob.30 <= .8), nudge_x = -.1, segment.alpha = .25) +
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
theta.9581 <- as.data.frame(cong.fit.b.rstan, c("theta_p[95]", "theta_p[81]"))
mean(theta.9581[, 1] > theta.9581[, 2])
# [1] 0.5447083
theta.8182 <- as.data.frame(cong.fit.b.rstan, c("theta_p[81]", "theta_p[82]"))
mean(theta.8182[, 1] > theta.8182[, 2])  # 81 only marginally higher probability than 82
# [1] 0.5157083

# APPENDIX CONTENT ----

# Rankplots for congeneric model, Gaussian
mcmc_rank_overlay(cong.fit.ln.rstan, regex_pars = c("alpha", "beta", "lambda", "sigma"),
                  ref_line = TRUE)
ggsave(paste0(print.images, "app_01_rank_gauss_cong.pdf"), height = 6, width = 8)
# Rankplots for congeneric model, Beta
mcmc_rank_overlay(cong.fit.b.rstan, regex_pars = c("i_means", "lambda", "sigma"),
                  ref_line = TRUE)
ggsave(paste0(print.images, "app_02_rank_beta_cong.pdf"), height = 6, width = 8)

# Posterior predictive for congeneric CFA as regression (Gaussian)
G.cong.ln <- as.data.frame(cong.fit.ln.rstan, paste0("yhat"))
# Should be 100 samples by 603 responses
dim(G.cong.ln <- G.cong.ln[sample(1:nrow(G.cong.ln), 100), ])

data.frame(item = X$item_n, response = X$response.01, t(G.cong.ln)) %>%
  pivot_longer(contains("X"), names_pattern = "X(.*)") %>%
  ggplot(aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(labels = percent_format(), name = "Distribution of prevalence") +
  geom_density(aes(value, group = name), fill = NA, col = "darkgrey") +
  geom_density(col = 1, fill = NA) +
  theme_classic() +
  theme(strip.background = element_blank())
ggsave(paste0(print.images, "app_03_ppc_gauss_cong.pdf"), height = 4, width = 6)

# More efficient estimation (of beta congeneric model) via "marginal likelihood" ----

# Compile Stan model
saveRDS(cmdstan_model(file.path(paste0(stan.scripts, "cong_ln_beta_sample_size_marg.stan"))),
        paste0(stan.scripts, "cong_ln_beta_sample_size_marg.rds"))
(cong.mod.b.m <- readRDS(paste0(stan.scripts, "cong_ln_beta_sample_size_marg.rds")))
# Same configuration settings as Gaussian model
cong.fit.b.m <- cong.mod.b.m$sample(
  data = data.b, seed = 12345, num_warmup = 1e3,
  num_samples = 1e3, num_chains = 4, num_cores = 4)

cong.fit.b.m$cmdstan_diagnose()
cong.fit.b.m.rstan <- read_stan_csv(cong.fit.b.m$output_files())

print(cong.fit.b.m.rstan, c("latent_mv"), include = FALSE)
#               mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# alpha        -0.92    0.01  0.40   -1.72   -1.17   -0.92   -0.67   -0.10  3149    1
# sigma_beta    0.88    0.01  0.32    0.47    0.67    0.81    1.01    1.69  2466    1
# beta[1]       0.02    0.00  0.19   -0.35   -0.10    0.02    0.15    0.40  4635    1
# beta[2]      -0.43    0.00  0.15   -0.73   -0.53   -0.43   -0.33   -0.13  4814    1
# beta[3]      -1.22    0.00  0.12   -1.47   -1.30   -1.22   -1.14   -0.99  3406    1
# beta[4]      -1.53    0.00  0.24   -2.00   -1.68   -1.52   -1.38   -1.06  5266    1
# beta[5]      -0.53    0.00  0.14   -0.82   -0.63   -0.53   -0.44   -0.26  4052    1
# beta[6]      -1.96    0.00  0.15   -2.25   -2.06   -1.96   -1.86   -1.66  4035    1
# lambda[1]     1.66    0.00  0.16    1.36    1.54    1.65    1.76    2.00  2300    1
# lambda[2]     1.34    0.00  0.12    1.11    1.25    1.33    1.41    1.59  4013    1
# lambda[3]     0.95    0.00  0.10    0.77    0.88    0.94    1.01    1.14  3198    1
# lambda[4]     2.24    0.00  0.17    1.92    2.12    2.23    2.35    2.60  4361    1
# lambda[5]     1.17    0.00  0.12    0.95    1.09    1.17    1.25    1.43  2763    1
# lambda[6]     1.23    0.00  0.11    1.02    1.15    1.22    1.30    1.47  3789    1
# prec[1]       5.73    0.02  0.92    4.09    5.09    5.65    6.31    7.70  2830    1
# prec[2]      10.80    0.03  1.67    7.75    9.61   10.71   11.88   14.28  4396    1
# prec[3]      13.39    0.03  1.99    9.67   11.98   13.35   14.68   17.51  3930    1
# prec[4]      69.13    0.66 21.13   35.96   54.05   66.48   81.41  116.92  1033    1
# prec[5]       6.85    0.02  0.99    5.05    6.17    6.80    7.47    8.96  3819    1
# prec[6]      17.80    0.05  2.96   12.49   15.73   17.61   19.64   24.14  3051    1
# i_means[1]    0.51    0.00  0.05    0.41    0.47    0.51    0.54    0.60  4629    1
# i_means[2]    0.39    0.00  0.04    0.32    0.37    0.39    0.42    0.47  4821    1
# i_means[3]    0.23    0.00  0.02    0.19    0.21    0.23    0.24    0.27  3396    1
# i_means[4]    0.18    0.00  0.03    0.12    0.16    0.18    0.20    0.26  5271    1
# i_means[5]    0.37    0.00  0.03    0.31    0.35    0.37    0.39    0.44  4065    1
# i_means[6]    0.12    0.00  0.02    0.10    0.11    0.12    0.14    0.16  3987    1
# sigma[1]      0.19    0.00  0.01    0.17    0.18    0.19    0.20    0.22  2599    1
# sigma[2]      0.14    0.00  0.01    0.12    0.14    0.14    0.15    0.17  4256    1
# sigma[3]      0.11    0.00  0.01    0.10    0.10    0.11    0.12    0.13  3369    1
# sigma[4]      0.05    0.00  0.01    0.03    0.04    0.05    0.05    0.07  1136    1
# sigma[5]      0.17    0.00  0.01    0.15    0.17    0.17    0.18    0.20  3652    1
# sigma[6]      0.08    0.00  0.01    0.06    0.07    0.08    0.08    0.10  2943    1
# lp__       1454.85    0.59 19.22 1416.57 1441.84 1455.39 1468.12 1490.92  1063    1

# R Session Info (Contains exact package versions and system information) ----
sessionInfo()  # Output is truncated to highlight important elements
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cmdstanr_0.0.0.9000     rstan_2.19.2            StanHeaders_2.19.0      bayesplot_1.7.0        
# [5] ggrepel_0.8.1           directlabels_2018.05.22 ggplot2_3.2.1           latex2exp_0.4.0        
# [9] patchwork_0.0.1         scales_1.1.1            tidyr_1.1.1             dplyr_1.0.2            
# [13] glmmTMB_1.0.0           lavaan_0.6-6            psych_2.0.7             data.table_1.12.6      
# 
# loaded via a namespace (and not attached):
#   [1] jsonlite_1.7.0     splines_3.6.0      assertthat_0.2.1   posterior_0.0.2    stats4_3.6.0       pbivnorm_0.6.0    
# [7] backports_1.1.8    pillar_1.4.6       lattice_0.20-41    glue_1.4.1         quadprog_1.5-8     digest_0.6.25     
# [13] checkmate_2.0.0    minqa_1.2.4        colorspace_1.4-1   sandwich_2.5-1     Matrix_1.2-18      plyr_1.8.6        
# [19] pkgconfig_2.0.3    purrr_0.3.3        xtable_1.8-4       mvtnorm_1.0-12     processx_3.4.2     lme4_1.1-21       
# [25] emmeans_1.5.0      tibble_3.0.1       farver_2.0.1       generics_0.0.2     ellipsis_0.3.0     TH.data_1.0-10    
# [31] withr_2.2.0        TMB_1.7.16         lazyeval_0.2.2     cli_2.0.2          mnormt_1.5-5       survival_3.2-3    
# [37] magrittr_1.5       crayon_1.3.4       estimability_1.3   ps_1.3.4           fansi_0.4.1        nlme_3.1-142      
# [43] MASS_7.3-52        pkgbuild_1.0.6     tools_3.6.0        loo_2.1.0          prettyunits_1.1.1  lifecycle_0.2.0   
# [49] matrixStats_0.55.0 multcomp_1.4-13    stringr_1.4.0      munsell_0.5.0      callr_3.3.2        packrat_0.5.0     
# [55] compiler_3.6.0     rlang_0.4.7        grid_3.6.0         nloptr_1.2.2.2     ggridges_0.5.1     rstudioapi_0.11   
# [61] labeling_0.3       boot_1.3-25        gtable_0.3.0       codetools_0.2-16   abind_1.4-5        inline_0.3.15     
# [67] reshape2_1.4.4     R6_2.4.1           gridExtra_2.3      zoo_1.8-8          stringi_1.4.6      parallel_3.6.0    
# [73] Rcpp_1.0.4.6       vctrs_0.3.2        tidyselect_1.1.0   coda_0.19-3       
