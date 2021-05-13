data.location <- "summary_text/mpi_dat/"  # location of csv dataset

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load packages (Check for packages that should be installed) ----
source("parallel_mod_mpi_cmd_load_packages.R")

# Pre-process data (Worth reviewing source file) ----
source("parallel_mod_mpi_cmd_prep_data.R")

stan.scripts <- "summary_text/stan_scripts/"  # location of stan scripts
print.images <- "summary_text/text_clean/"  # location to print images

# FIGURE 1 ----
ggplot(X, aes(response.l)) + geom_histogram(col = 1) + facet_wrap(~ item_n) +
  scale_x_continuous(name = "Probit-transformed prevalence of deprivation") +
  theme(strip.background = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(.5, "cm"), panel.grid.minor.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(y = "Number of countries")
ggsave(paste0(print.images, "01_deprivation.pdf"), width = 6, height = 2.75)

# Table 1 ----
describe(dat[, 15:20])
#    vars   n  mean   sd median trimmed  mad   min  max range  skew kurtosis   se
# x1    1 100 -0.06 1.45  -0.17   -0.11 1.64 -2.34 2.58  4.91  0.17    -1.01 0.14
# x2    2 101  0.31 1.03   0.21    0.30 1.20 -1.95 2.24  4.19  0.10    -1.01 0.10
# x3    3 101  0.80 0.77   0.67    0.77 0.78 -0.53 2.72  3.26  0.36    -0.83 0.08
# x4    4 100  0.85 1.29   0.86    0.90 1.86 -1.66 2.58  4.24 -0.16    -1.37 0.13
# x5    5 100  0.38 1.04   0.33    0.40 1.21 -2.10 2.40  4.50 -0.19    -0.65 0.10
# x6    6 101  1.31 0.96   1.13    1.29 1.09 -0.53 3.50  4.02  0.17    -0.89 0.10

round(cor(dat[, 15:20], use = "p"), 2)

# MODELLING ----

# CONGENERIC CFA ----
# Congeneric CFA in lavaan ----
# Creating error SD as square root of error variance
X.tmp <- dat[, 15:20]
colnames(X.tmp) <- c("cf", "sanit", "dw", "elect", "hous", "ast")
writeLines(lav.form <- paste(
  "F =~ cf + sanit + dw + elect + hous + ast
  cf ~~ a * cf
  sanit ~~ b * sanit
  dw ~~ c * dw
  elect ~~ d * elect
  hous ~~ e * hous
  ast ~~ f * ast",
  "a1 := sqrt(a)\nb1 := sqrt(b)\nc1 := sqrt(c)\nd1 := sqrt(d)\ne1 := sqrt(e)\nf1 := sqrt(f)",
  sep = "\n"))
summary(lav.fit <- cfa(
  lav.form, X.tmp, std.lv = TRUE, meanstructure = TRUE, missing = "ML"),
  fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)

# FIGURE 2 ----
# plot fit
png(paste0(print.images, "02_lavfit.png"), width = 5, height = 2, units = "in", res = 150)
semPaths(lav.fit, "est", fade = FALSE, style = "ram", edge.label.cex = 1.5,
         mar = c(2, 2, 2, 2), edge.color = cbPalette[1])
dev.off()

parameterEstimates(lav.fit, standardized = TRUE, rsquare = TRUE)
#      lhs op     rhs label    est    se       z pvalue ci.lower ci.upper std.lv std.all std.nox
# 1      F =~      cf        1.280 0.113  11.361  0.000    1.059    1.501  1.280   0.888   0.888
# 2      F =~   sanit        0.926 0.078  11.791  0.000    0.772    1.079  0.926   0.906   0.906
# 3      F =~      dw        0.623 0.063   9.897  0.000    0.499    0.746  0.623   0.812   0.812
# 4      F =~   elect        1.248 0.093  13.453  0.000    1.066    1.429  1.248   0.974   0.974
# 5      F =~    hous        0.883 0.083  10.684  0.000    0.721    1.044  0.883   0.855   0.855
# 6      F =~     ast        0.882 0.073  12.142  0.000    0.740    1.025  0.882   0.921   0.921
# 14    cf ~1               -0.043 0.144 -0.296  0.767   -0.324    0.239 -0.043  -0.029  -0.029
# 15 sanit ~1                0.313 0.102  3.078  0.002    0.114    0.512  0.313   0.306   0.306
# 16    dw ~1                0.800 0.076 10.490  0.000    0.651    0.949  0.800   1.044   1.044
# 17 elect ~1                0.849 0.127  6.661  0.000    0.599    1.099  0.849   0.663   0.663
# 18  hous ~1                0.388 0.103  3.767  0.000    0.186    0.589  0.388   0.375   0.375
# 19   ast ~1                1.307 0.095 13.710  0.000    1.120    1.493  1.307   1.364   1.364
# 21    a1 := sqrt(a)    a1  0.664 0.052  12.700  0.000    0.562    0.766  0.664   0.460   0.460
# 22    b1 := sqrt(b)    b1  0.433 0.035  12.314  0.000    0.364    0.502  0.433   0.424   0.424
# 23    c1 := sqrt(c)    c1  0.447 0.033  13.428  0.000    0.382    0.512  0.447   0.583   0.583
# 24    d1 := sqrt(d)    d1  0.288 0.043   6.753  0.000    0.204    0.372  0.288   0.225   0.225
# 25    e1 := sqrt(e)    e1  0.536 0.041  13.048  0.000    0.456    0.617  0.536   0.519   0.519
# 26    f1 := sqrt(f)    f1  0.372 0.031  11.851  0.000    0.311    0.434  0.372   0.389   0.389
# 27    cf r2      cf        0.788    NA      NA     NA       NA       NA     NA      NA      NA
# 28 sanit r2   sanit        0.821    NA      NA     NA       NA       NA     NA      NA      NA
# 29    dw r2      dw        0.660    NA      NA     NA       NA       NA     NA      NA      NA
# 30 elect r2   elect        0.949    NA      NA     NA       NA       NA     NA      NA      NA
# 31  hous r2    hous        0.730    NA      NA     NA       NA       NA     NA      NA      NA
# 32   ast r2     ast        0.849    NA      NA     NA       NA       NA     NA      NA      NA

rm(X.tmp)

# FIGURE 3 ----
# plot log-normal
(ggplot(data.frame(x = seq(0, 2, 1e-4), density = dlnorm(seq(0, 2, 1e-4), log(.5), log(sqrt(3) / .5) / qnorm(.95))),
        aes(x, density)) + geom_line() + theme_classic() +
   scale_x_continuous(breaks = c(0, .5, sqrt(3)), labels = number_format()) +
   geom_vline(xintercept = c(.5, sqrt(3)), linetype = 2, alpha = .5) +
   theme(axis.title = element_blank()) +
   labs(tag = "A", subtitle = "Density")) +
  (ggplot(data.frame(x = seq(0, 2, 1e-4), cum_dens = plnorm(seq(0, 2, 1e-4), log(.5), log(sqrt(3) / .5) / qnorm(.95))),
          aes(x, cum_dens)) + geom_line() + theme_classic() +
     scale_x_continuous(breaks = c(0, .5, sqrt(3)), labels = number_format()) +
     scale_y_continuous(labels = percent_format(), breaks = c(.0, .25, .5, .75, .95)) +
     theme(axis.title.x = element_text(hjust = 1.5), axis.title.y = element_blank()) +
     geom_hline(yintercept = c(.5, .95), linetype = 2, alpha = .5) +
     geom_vline(xintercept = c(.5, sqrt(3)), linetype = 2, alpha = .5) +
     labs(x = "A-priori expectations for lambda_i parameter",
          subtitle = "Cumulative distribution", tag = "B"))
ggsave(paste0(print.images, "03_log_normal_lambda.pdf"), width = 5, height = 2.5)

# Congeneric regression in Stan ----
# Data list with references to equation 11
data <- list(
  alpha_scale = 3,  # scale of a parameter
  beta_scale = 2,  # scale of sigma_t
  sigma_scale = sqrt(3),  # scale of sigma_i
  lambda_median = log(.5),  # log-median of lambda_i
  lambda_scale = log(sqrt(3) / .5) / qnorm(.95),  # scale of lambda_i
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

# Compile Stan model (takes time)
mod <- cmdstan_model(file.path(paste0(stan.scripts, "cong.stan")))
fit.0 <- mod$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems,
)
fit.0$cmdstan_diagnose()  # Run quick diagnostic checks
fit.0 <- read_stan_csv(fit.0$output_files())  # convert to Rstan

# Print major parameters
print(fit.0, c("theta_p", "yhat", "lp__"), include = FALSE)
#              mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha        0.56    0.01 0.79 -1.01  0.03  0.56  1.09  2.12  2994 1.00
# beta_dev[1] -0.59    0.01 0.79 -2.16 -1.12 -0.60 -0.06  0.96  3101 1.00
# beta_dev[2] -0.24    0.01 0.79 -1.80 -0.77 -0.24  0.29  1.31  3098 1.00
# beta_dev[3]  0.25    0.01 0.79 -1.31 -0.28  0.24  0.78  1.80  3077 1.00
# beta_dev[4]  0.30    0.01 0.79 -1.26 -0.23  0.29  0.82  1.83  3097 1.00
# beta_dev[5] -0.16    0.01 0.79 -1.73 -0.69 -0.17  0.36  1.38  3104 1.00
# beta_dev[6]  0.75    0.01 0.79 -0.81  0.22  0.75  1.28  2.30  3088 1.00
# lambda[1]    1.24    0.00 0.11  1.04  1.17  1.24  1.31  1.47  1943 1.00
# lambda[2]    0.90    0.00 0.07  0.76  0.85  0.90  0.95  1.06  1849 1.00
# lambda[3]    0.60    0.00 0.06  0.49  0.56  0.60  0.64  0.73  2461 1.00
# lambda[4]    1.22    0.00 0.09  1.06  1.16  1.21  1.27  1.40  1331 1.01
# lambda[5]    0.86    0.00 0.08  0.71  0.80  0.85  0.91  1.02  2041 1.00
# lambda[6]    0.86    0.00 0.07  0.73  0.81  0.86  0.90  1.01  1669 1.00
# sigma[1]     0.68    0.00 0.05  0.58  0.64  0.68  0.72  0.80 17469 1.00
# sigma[2]     0.44    0.00 0.04  0.38  0.42  0.44  0.47  0.52 16042 1.00
# sigma[3]     0.46    0.00 0.03  0.40  0.43  0.46  0.48  0.53 19362 1.00
# sigma[4]     0.29    0.00 0.05  0.20  0.26  0.29  0.32  0.38  3763 1.00
# sigma[5]     0.55    0.00 0.04  0.47  0.52  0.55  0.58  0.64 17561 1.00
# sigma[6]     0.38    0.00 0.03  0.32  0.36  0.38  0.40  0.45 13651 1.00
# beta[1]     -0.04    0.00 0.14 -0.31 -0.13 -0.04  0.06  0.24  1140 1.00
# beta[2]      0.32    0.00 0.10  0.13  0.25  0.32  0.38  0.51  1131 1.00
# beta[3]      0.80    0.00 0.07  0.66  0.75  0.80  0.85  0.95  1356 1.00
# beta[4]      0.86    0.00 0.12  0.61  0.77  0.85  0.94  1.10   891 1.01
# beta[5]      0.39    0.00 0.10  0.20  0.32  0.39  0.46  0.59  1249 1.00
# beta[6]      1.31    0.00 0.09  1.13  1.25  1.31  1.37  1.49  1027 1.00

# saveRDS(fit.0, "res/mod_1.rds")
fit.0 <- readRDS("res/mod_1.rds")
check_hmc_diagnostics(fit.0)

# coefficient omega:
sum(apply(as.data.frame(fit.0, "lambda"), 2, median)) ^ 2 / (
  sum(apply(as.data.frame(fit.0, "lambda"), 2, median)) ^ 2 +
    sum(apply(as.data.frame(fit.0, "sigma"), 2, median) ^ 2)
)

# Traceplots:
mcmc_trace(fit.0, pars = c(
  paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
  paste0("sigma[", 1:6, "]")))
# Rank plots:
mcmc_rank_overlay(
  fit.0, pars = c(
    paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
    paste0("sigma[", 1:6, "]")), ref_line = TRUE) +
  scale_x_continuous(labels = number_format(scale = 1e-3, 1), name = "Rank (1000s)") +
  theme(legend.position = "top", strip.background = element_blank())
ggsave(paste0(print.images, "app_01_rank_mod_1.pdf"), width = 6.5, height = 5)

pars.0.df <- as.data.frame(
  summary(fit.0, c("beta", "lambda", "sigma"))$summary)
# beta[1-6] are means
# lambda[1-6] are loadings
# sigma[1-6] are residual standard deviations
pars.0.df$variable <- rownames(pars.0.df)
# Join lavaan and Bayes results
pars.0.df
pars.0.df.joined <- as.data.table(cbind(
  pars.0.df[, c(11, 6, 4, 8)],
  parameterEstimates(lav.fit)[c(14:19, 1:6, 21:26), c(5, 9:10)]))
pars.0.df.joined
# Rename columns
setnames(pars.0.df.joined, 2:7, c("m.1", "ll.1", "ul.1", "m.2", "ll.2", "ul.2"))
pars.0.df.joined <- reshape(
  pars.0.df.joined, direction = "long", varying = list(c(2, 5), c(3, 6), c(4, 7)))
pars.0.df.joined[, model := c("Bayesian", "ML")[time]]
pars.0.df.joined[, type := factor(
  gsub("\\[\\d\\]", "", variable), levels = c("lambda", "beta", "sigma"),
  labels = c("Loadings", "Means", "Residual SD"))]
pars.0.df.joined[, item := as.integer(gsub("[a-z]|\\[|\\]", "", variable))]
pars.0.df.joined

fig.capt <-
  "Indicator codes: 1 = Cooking fuel, 2 = Sanitation, 3 = Drinking water, 4 = Electricity, 5 = Housing, 6 = Assets"
ggplot(pars.0.df.joined, aes(item, m.1, fill = model, shape = model)) +
  geom_pointrange(aes(ymin = ll.1, ymax = ul.1), position = position_dodge(.5), alpha = .65) +
  theme_classic() + theme(legend.position = "top", strip.background = element_blank(),
                          axis.title.y = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ type, scales = "free_y") +
  labs(x = "Indicator number", subtitle = "Estimates with 95% interval",
       shape = "Method", fill = "Method", caption = fig.capt)
ggsave(paste0(print.images, "04_params_lav_bayes.pdf"), width = 6.5, height = 4)

rm(pars.0.df, pars.0.df.joined)
gc()

# Extract predicted outcomes
dim(yhat.0.f <- as.data.frame(fit.0, "yhat"))

yhat.0.df <- as.data.table(t(apply(yhat.0.f, 2, quantile, c(.5, .25, .75, .025, .975))))
yhat.0.df$item <- X$item
yhat.0.df$item_n <- X$item_n
yhat.0.df$id <- X$id
yhat.0.df$response <- X$response.l

setnames(yhat.0.df, c("50%", "25%", "75%", "2.5%", "97.5%"), c("median", "lli", "uli", "ll", "ul"))
yhat.0.df

yhat.0.df[, mean(response > ll & response < ul), item_n]
#            item_n        V1
# 1:   Cooking fuel 0.9700000
# 2:     Sanitation 0.9306931
# 3: Drinking water 0.9504950
# 4:    Electricity 1.0000000
# 5:        Housing 0.9500000
# 6:         Assets 0.9603960

yhat.0.df$theta <- apply(as.data.frame(fit.0, "theta_p"), 2, median)[yhat.0.df$id]

# Heteroskedastic pattern present (for cooking fuel and drinking water)
ggplot(yhat.0.df, aes(theta, response - median)) +
  geom_point(shape = 1, alpha = .25) +
  geom_abline(linetype = 2, intercept = 0, slope = 0) +
  facet_wrap(~ item_n) +
  geom_smooth(formula = y ~ s(x, bs = "ts"), method = "gam", col = cbbPalette[1],
              se = TRUE, size = .5) +
  geom_smooth(aes(y = response - ll), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  geom_smooth(aes(y = response - ul), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  labs(x = "Standardized country average (k_c^')", y = "residual with smoothed 95% interval", col = "Model")
ggsave(paste0(print.images, "05_resid_bayes.pdf"), width = 6.5, height = 4.5)

# Should be 100 samples by 603 responses
set.seed(12345)
dim(yhat.0 <- yhat.0.f[sample(1:nrow(yhat.0.f), 100), ])
yhat.0 <- data.table(item = X$item_n, response = X$response.l, t(yhat.0))
setnames(yhat.0, 3:ncol(yhat.0), paste0("post.", 1:(ncol(yhat.0) - 2)))
yhat.0 <- reshape(yhat.0, direction = "long", varying = 3:ncol(yhat.0))
yhat.0

ggplot(yhat.0, aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "#999999") +
  geom_density(col = 1, fill = NA, data = yhat.0[time == 1]) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

# Compile Stan model (takes time)
mod.1 <- cmdstan_model(file.path(paste0(stan.scripts, "cong_joint.stan")))
fit.1 <- mod.1$sample(
  data = data,  # data list passed to Stan
  seed = 12345,  # Seed for reproducibility
  iter_warmup = 1e3,  # 1000 warmup samples per chain
  iter_sampling = 3e3,  # 3000 samples for inference per chain
  chains = 8,  # 8 chains
  parallel_chains = 8  # 8 cores on multicore systems
)
fit.1$cmdstan_diagnose()  # Run quick diagnostic checks
fit.1 <- read_stan_csv(fit.1$output_files())  # convert to Rstan

# Print major parameters
print(fit.1, c("Theta", "yhat", "lp__"), include = FALSE)
#                   mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha             0.57    0.01 0.79 -0.99  0.03  0.57  1.10  2.13  4625 1.00
# beta_dev[1]      -0.59    0.01 0.79 -2.15 -1.12 -0.60 -0.06  0.94  4683 1.00
# beta_dev[2]      -0.25    0.01 0.79 -1.80 -0.78 -0.26  0.28  1.29  4712 1.00
# beta_dev[3]       0.21    0.01 0.79 -1.34 -0.32  0.20  0.74  1.75  4716 1.00
# beta_dev[4]       0.28    0.01 0.79 -1.27 -0.25  0.27  0.81  1.81  4693 1.00
# beta_dev[5]      -0.18    0.01 0.79 -1.74 -0.71 -0.19  0.35  1.36  4725 1.00
# beta_dev[6]       0.73    0.01 0.79 -0.82  0.21  0.73  1.26  2.27  4714 1.00
# lambda[1]         1.28    0.00 0.11  1.08  1.20  1.27  1.35  1.50  2565 1.00
# lambda[2]         0.93    0.00 0.08  0.79  0.88  0.93  0.98  1.09  2685 1.00
# lambda[3]         0.62    0.00 0.06  0.51  0.58  0.61  0.66  0.74  3744 1.00
# lambda[4]         1.21    0.00 0.09  1.05  1.15  1.21  1.27  1.39  2130 1.00
# lambda[5]         0.86    0.00 0.08  0.72  0.81  0.86  0.91  1.02  3255 1.00
# lambda[6]         0.86    0.00 0.07  0.73  0.81  0.85  0.90  1.00  2560 1.00
# ln_alpha         -1.73    0.01 0.29 -2.30 -1.88 -1.73 -1.57 -1.18  2837 1.00
# ln_sigma2_dev[1]  0.62    0.01 0.34  0.01  0.40  0.60  0.82  1.32  3397 1.00
# ln_sigma2_dev[2] -0.07    0.01 0.33 -0.72 -0.26 -0.07  0.13  0.59  3296 1.00
# ln_sigma2_dev[3] -0.14    0.01 0.33 -0.78 -0.33 -0.14  0.06  0.49  3333 1.00
# ln_sigma2_dev[4] -0.58    0.01 0.36 -1.34 -0.79 -0.55 -0.34  0.08  4287 1.00
# ln_sigma2_dev[5]  0.40    0.01 0.32 -0.20  0.20  0.38  0.58  1.06  3355 1.00
# ln_sigma2_dev[6] -0.21    0.01 0.32 -0.83 -0.40 -0.21 -0.03  0.41  3397 1.00
# ls_scale          0.59    0.00 0.26  0.24  0.41  0.54  0.71  1.25  5766 1.00
# ll_scale          0.64    0.00 0.25  0.29  0.47  0.60  0.76  1.25 10816 1.00
# ln_lambda_pos[1]  0.69    0.00 0.24  0.27  0.52  0.67  0.84  1.20  5193 1.00
# ln_lambda_pos[2]  0.76    0.00 0.22  0.32  0.60  0.75  0.90  1.22  9892 1.00
# ln_lambda_oth[1] -0.48    0.01 0.32 -1.16 -0.68 -0.46 -0.25  0.10  3278 1.00
# ln_lambda_oth[2]  0.00    0.00 0.27 -0.53 -0.17  0.00  0.18  0.53  9098 1.00
# ln_lambda_oth[3] -0.33    0.00 0.27 -0.84 -0.50 -0.33 -0.16  0.21  7216 1.00
# ln_lambda_oth[4]  0.12    0.00 0.27 -0.39 -0.07  0.11  0.29  0.68  6218 1.00
# r_tr              0.82    0.00 0.07  0.67  0.77  0.82  0.87  0.96  1614 1.01
# ln_lambda[1]      0.69    0.00 0.24  0.27  0.52  0.67  0.84  1.20  5193 1.00
# ln_lambda[2]     -0.48    0.01 0.32 -1.16 -0.68 -0.46 -0.25  0.10  3278 1.00
# ln_lambda[3]      0.76    0.00 0.22  0.32  0.60  0.75  0.90  1.22  9892 1.00
# ln_lambda[4]      0.00    0.00 0.27 -0.53 -0.17  0.00  0.18  0.53  9098 1.00
# ln_lambda[5]     -0.33    0.00 0.27 -0.84 -0.50 -0.33 -0.16  0.21  7216 1.00
# ln_lambda[6]      0.12    0.00 0.27 -0.39 -0.07  0.11  0.29  0.68  6218 1.00
# beta[1]          -0.03    0.00 0.14 -0.30 -0.12 -0.03  0.07  0.25  1768 1.01
# beta[2]           0.32    0.00 0.10  0.12  0.25  0.32  0.38  0.52  1714 1.01
# beta[3]           0.78    0.00 0.08  0.63  0.72  0.77  0.83  0.93  2120 1.01
# beta[4]           0.84    0.00 0.12  0.60  0.76  0.84  0.93  1.09  1564 1.01
# beta[5]           0.38    0.00 0.10  0.18  0.32  0.38  0.45  0.58  1996 1.01
# beta[6]           1.30    0.00 0.09  1.12  1.24  1.30  1.36  1.48  1715 1.01
# sigma[1]          0.58    0.00 0.06  0.46  0.54  0.58  0.62  0.70  6870 1.00
# sigma[2]          0.41    0.00 0.04  0.33  0.38  0.41  0.44  0.49  6281 1.00
# sigma[3]          0.39    0.00 0.04  0.32  0.37  0.39  0.42  0.47  9874 1.00
# sigma[4]          0.32    0.00 0.04  0.24  0.29  0.32  0.35  0.40  9416 1.00
# sigma[5]          0.52    0.00 0.05  0.43  0.48  0.51  0.54  0.61 12085 1.00
# sigma[6]          0.38    0.00 0.03  0.32  0.36  0.38  0.40  0.45 21225 1.00
# r                 0.64    0.00 0.15  0.33  0.55  0.65  0.75  0.91  1614 1.01

# print(describe(apply(as.data.frame(fit.0, "theta_p"), 2, mad)), digits = 4)
# print(describe(apply(as.data.frame(fit.1, "Theta"), 2, mad)[1:101]), digits = 4)
# hist(apply(as.data.frame(fit.0, "theta_p"), 2, mad))
# hist(apply(as.data.frame(fit.1, "Theta"), 2, mad)[1:101])

# saveRDS(fit.1, "res/mod_2.rds")
fit.1 <- readRDS("res/mod_2.rds")
check_hmc_diagnostics(fit.1)
# 1 of 24000 iterations ended with a divergence (0.00416666666666667%).

# coefficient omega:
sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 / (
  sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 +
    sum(apply(as.data.frame(fit.1, "sigma"), 2, median) ^ 2)
)

apply(as.data.frame(fit.1, "lambda"), 2, median)
apply(as.data.frame(fit.1, "sigma"), 2, median)

hist(sapply(1:101, function (i) {
  l.sum.2 <- sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2
  l.sum.2 / (
    l.sum.2 +
      sum(exp((apply(as.data.frame(fit.1, "ln_alpha"), 2, median) +
                 apply(as.data.frame(fit.1, "ln_sigma2_dev"), 2, median) +
                 apply(as.data.frame(fit.1, "ln_lambda"), 2, median) *
                 apply(as.data.frame(fit.1, paste0("Theta[", i, ",2]")), 2, median))))
  )
}))
hist(sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 / (
  sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 +
    aggregate(apply(as.data.frame(fit.1, "yhat"), 2, var) ~ X$id, FUN = sum)[, 2]
))
sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 / (
  sum(apply(as.data.frame(fit.1, "lambda"), 2, median)) ^ 2 +
    sum(exp((apply(as.data.frame(fit.1, "ln_alpha"), 2, median) +
               apply(as.data.frame(fit.1, "ln_sigma2_dev"), 2, median) +
               apply(as.data.frame(fit.1, "ln_lambda"), 2, median) *
               apply(as.data.frame(fit.1, "Theta"), 2, median)[103])))
)

# Traceplots:
mcmc_trace(fit.1, pars = c(
  paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
  paste0("sigma[", 1:6, "]"), paste0("ln_lambda[", 1:6, "]")))
# Rank plots:
mcmc_rank_overlay(
  fit.1, pars = c(
    paste0("beta[", 1:6, "]"), paste0("lambda[", 1:6, "]"),
    paste0("sigma[", 1:6, "]"), paste0("ln_lambda[", 1:6, "]")),
  ref_line = TRUE)

{
  pars.1.df <- as.data.frame(summary(fit.1, c("ln_lambda"))$summary)
  pars.1.df$variable <- rownames(pars.1.df)
  pars.1.df
  pars.1.df <- as.data.table(pars.1.df[, c(4, 6, 8, 11)])
  pars.1.df
  setnames(pars.1.df, 1:3, c("ll", "median", "ul"))
  pars.1.df[, item := as.integer(gsub("[a-z]|\\[|\\]|_", "", variable))]
  pars.1.df
  pars.1.df$item_t <- levels(X$item_n)
  pars.1.df
}

plt.pars <- ggplot(pars.1.df, aes(reorder(item_t, -item), median)) + coord_flip() +
  geom_pointrange(aes(ymin = ll, ymax = ul), position = position_dodge(.5), alpha = .65, shape = 4) +
  theme_bw() + theme(legend.position = "top", strip.background = element_blank(),
                     axis.title.y = element_blank(), panel.border = element_blank(),
                     axis.ticks = element_blank()) +
  scale_shape_manual(values = c(1, 4)) +
  labs(x = "Indicator number", y = "Log-variance loadings with 95% quantile interval")
plt.pars

rm(pars.1.df)
gc()

dim(yhat.1.f <- as.data.frame(fit.1, "yhat"))

yhat.1.df <- as.data.table(t(apply(yhat.1.f, 2, quantile, c(.5, .25, .75, .025, .975))))
yhat.1.df$item <- X$item
yhat.1.df$item_n <- X$item_n
yhat.1.df$id <- X$id
yhat.1.df$response <- X$response.l

setnames(yhat.1.df, c("50%", "25%", "75%", "2.5%", "97.5%"), c("median", "lli", "uli", "ll", "ul"))
yhat.1.df

yhat.1.df[, mean(response > ll & response < ul), item_n]
#            item_n       V1
# 1:   Cooking fuel 0.980000
# 2:     Sanitation 0.990099
# 3: Drinking water 0.980198
# 4:    Electricity 1.000000
# 5:        Housing 0.950000
# 6:         Assets 0.960396

yhat.1.df$theta <- colMeans(as.data.frame(fit.1, "Theta"))[1:101][yhat.1.df$id]

# Heteroskedastic pattern present (for cooking fuel and drinking water)
plt.resid <- ggplot(yhat.1.df, aes(theta, response - median)) +
  geom_point(shape = 1, alpha = .25) +
  geom_abline(linetype = 2, intercept = 0, slope = 0) +
  facet_wrap(~ item_n, scales = "free_x") +
  geom_smooth(formula = y ~ s(x, bs = "ts"), method = "gam", col = cbbPalette[1],
              se = TRUE, size = .5) +
  geom_smooth(aes(y = response - ll), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  geom_smooth(aes(y = response - ul), formula = y ~ s(x, bs = "ts"),
              method = "gam", col = cbbPalette[1], se = FALSE, size = .5) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  labs(x = "Country averages", y = "residual with smoothed 95% interval", col = "Model")
(plt.pars + labs(tag = "A")) / (plt.resid + labs(tag = "B")) + plot_layout(heights = c(1.5, 4))
ggsave(paste0(print.images, "06_resid_bayes.pdf"), width = 6.5, height = 5.5)

# Should be 100 samples by 603 responses
set.seed(12345)
dim(yhat.1 <- yhat.1.f[sample(1:nrow(yhat.1.f), 100), ])
yhat.1 <- data.table(item = X$item_n, response = X$response.l, t(yhat.1))
setnames(yhat.1, 3:ncol(yhat.1), paste0("post.", 1:(ncol(yhat.1) - 2)))
yhat.1 <- reshape(yhat.1, direction = "long", varying = 3:ncol(yhat.1))
yhat.1

ggplot(yhat.1, aes(response)) + facet_wrap(~ item, scales = "free") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "darkgrey") +
  geom_density(aes(response), col = 1, fill = NA, data = yhat.1[time == 1]) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

ppc.joined <- rbindlist(list(yhat.0, yhat.1), idcol = "model")
ppc.joined[, model.t := c("Location", "Location-scale")[model]]

ggplot(ppc.joined, aes(response)) +
  facet_grid_paginate(item ~ model.t, scales = "free_y") +
  scale_x_continuous(name = "Distribution of data against posterior draws") +
  geom_density(aes(post, group = factor(time)), fill = NA, col = "darkgrey", size = .25) +
  geom_density(col = 1, fill = NA, data = yhat.1[time == 1]) +
  coord_cartesian(xlim = c(-10, 10)) +
  theme_classic() +
  theme(strip.background = element_blank(), plot.margin = margin(0, 0, 0, 0))

rm(yhat.0, yhat.1, ppc.joined)
gc()

yhat.df <- rbindlist(list(yhat.0.df, yhat.1.df), idcol = "model")
yhat.df
rm(yhat.0.df, yhat.1.df)
gc()

yhat.df[, model.t := c("Basic", "Location-scale")[model]]

Theta <- rbindlist(list(
  data.table(
    id = as.integer(gsub("_|[a-z]+|\\[|\\]", "", colnames(as.data.frame(fit.0, "theta_p")))),
    theta = colMeans(as.data.frame(fit.0, "theta_p")), model = 1),
  data.table(
    id = as.integer(gsub("T|,|[a-z]+|\\[|1\\]", "", colnames(as.data.frame(fit.1, "Theta"))[1:101])),
    theta = colMeans(as.data.frame(fit.1, "Theta"))[1:101], model = 2)))

yhat.df <- merge(yhat.df, Theta)
yhat.df[, theta.m := mean(theta), id]
yhat.df

ggplot(yhat.df[model == 2], aes(theta, median)) +
  geom_point(aes(theta, response), shape = 1, alpha = .125, size = .5) +
  geom_ribbon(aes(ymin = ll, ymax = ul), alpha = .125, col = cbPalette[1]) +
  geom_line(alpha = .5, linetype = 2) +
  geom_smooth(aes(y = ll), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .25, col = cbbPalette[1], linetype = 2) +
  geom_smooth(aes(y = ul), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .25, col = cbbPalette[1], linetype = 2) +
  facet_wrap(item_n ~ ., scales = "fixed") +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  scale_color_manual(values = cbbPalette[-1]) +
  labs(x = "Standardized country averages", y = "95% quantile interval",
       col = "Model", linetype = "Model")
ggsave(paste0(print.images, "07_model_suggestions_1.pdf"), width = 6.5, height = 4)

ggplot(yhat.df, aes(theta, median)) +
  geom_point(aes(theta, response), shape = 1, alpha = .25, size = .5) +
  geom_line(alpha = .5) +
  geom_smooth(aes(y = ll), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .5, col = cbbPalette[1]) +
  geom_smooth(aes(y = ul), formula = y ~ s(x, bs = "ts"), method = "gam",
              se = FALSE, size = .5, col = cbbPalette[1]) +
  facet_wrap(item_n ~ ., scales = "fixed", nrow = 1) +
  theme_classic() + theme(strip.background = element_blank(), legend.position = "top") +
  scale_color_manual(values = cbbPalette[-1]) +
  labs(x = "Standardized country averages", y = "Smoothed 95% quantile interval",
       col = "Model", linetype = "Model")
ggsave(paste0(print.images, "07_model_suggestions.pdf"), width = 6.5, height = 2.5)

# R Session Info (Contains exact package versions and system information) ----
sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 20.04 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
# LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] cmdstanr_0.2.2         rstan_2.21.3           StanHeaders_2.21.0-7   bayesplot_1.7.2        ggrepel_0.8.2         
# [6] directlabels_2020.6.17 ggforce_0.3.2          ggplot2_3.3.3          latex2exp_0.4.0        patchwork_1.1.0       
# [11] scales_1.1.1           lavaan_0.6-7           psych_2.0.9            data.table_1.14.0     
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6         lattice_0.20-41    prettyunits_1.1.1  ps_1.5.0           digest_0.6.27      assertthat_0.2.1  
# [7] utf8_1.2.1         V8_3.4.0           R6_2.5.0           plyr_1.8.6         backports_1.2.1    ggridges_0.5.2    
# [13] stats4_4.0.3       pillar_1.6.0       rlang_0.4.11       curl_4.3.1         callr_3.5.1        Matrix_1.2-18     
# [19] checkmate_2.0.0    pbivnorm_0.6.0     splines_4.0.3      labeling_0.4.2     stringr_1.4.0      loo_2.4.1         
# [25] polyclip_1.10-0    munsell_0.5.0      compiler_4.0.3     pkgconfig_2.0.3    mnormt_2.0.2       pkgbuild_1.2.0    
# [31] tmvnsim_1.0-2      mgcv_1.8-33        tidyselect_1.1.1   tibble_3.1.1       gridExtra_2.3      codetools_0.2-16  
# [37] matrixStats_0.58.0 quadprog_1.5-8     fansi_0.4.2        crayon_1.4.1       dplyr_1.0.5        withr_2.4.2       
# [43] MASS_7.3-53        grid_4.0.3         nlme_3.1-149       jsonlite_1.7.2     gtable_0.3.0       lifecycle_1.0.0   
# [49] DBI_1.1.0          magrittr_2.0.1     RcppParallel_5.0.2 cli_2.5.0          stringi_1.5.3      farver_2.1.0      
# [55] ellipsis_0.3.2     generics_0.1.0     vctrs_0.3.8        tools_4.0.3        glue_1.4.2         tweenr_1.0.1      
# [61] purrr_0.3.4        processx_3.4.5     parallel_4.0.3     inline_0.3.17      colorspace_2.0-1  

