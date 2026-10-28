# Run from the root of the flexbayes package.
# This script generates datasets only; no Bayesian model fitting.
set.seed(2323)
n = 80
x = rnorm(n)
z = 0.4 * x + rnorm(n, sd = sqrt(1 - 0.4^2))
a = factor(sample(c("control", "treatment"), n, replace = TRUE))
b = factor(sample(c("low", "medium", "high"), n, replace = TRUE))
y = 2 + 0.7 * x - 0.3 * z + 0.5 * (a == "treatment") +
  0.4 * x * (a == "treatment") + rnorm(n, sd = 1)
prob = plogis(-0.4 + 0.9 * x - 0.6 * z + 0.5 * (a == "treatment"))
y_bin = rbinom(n, size = 1, prob = prob)
y_pois = rpois(n, lambda = exp(0.3 + 0.35 * x - 0.2 * z))
small_bayes = data.frame(y, x, z, a, b, y_bin, y_pois)
usethis::use_data(small_bayes, overwrite = TRUE)

# Business example: threshold for the advertising slope is 1.
# The fitted probability of exceeding 1 depends on sampling and priors.
set.seed(2324)
n = 90
advertising = runif(n, min = 5, max = 45)
revenue = 60 + 1.13 * advertising + rnorm(n, sd = 17)
marketing_bayes = data.frame(market_id = seq_len(n), advertising, revenue)
usethis::use_data(marketing_bayes, overwrite = TRUE)

# Grouped data for future hierarchical model testing.
set.seed(2325)
n_groups = 8
n_each = 12
id = factor(rep(paste0("site", seq_len(n_groups)), each = n_each))
x = rnorm(length(id))
a = factor(sample(c("control", "treatment"), length(id), replace = TRUE))
site_effect = rnorm(n_groups, sd = 0.7)
y = 1 + 0.6 * x + 0.5 * (a == "treatment") +
  site_effect[as.integer(id)] + rnorm(length(id), sd = 0.8)
small_bayes_mixed = data.frame(y, x, a, id)
usethis::use_data(small_bayes_mixed, overwrite = TRUE)
