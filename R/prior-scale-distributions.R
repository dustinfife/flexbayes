# Positive, zero-centered halves match rstanarm's Gaussian prior_aux families.
prior_exp_survival = function(q, rate) stats::pexp(q, rate, lower.tail = FALSE)
prior_dhn = function(x, scale) ifelse(x < 0, 0, 2 * stats::dnorm(x, 0, scale))
prior_shn = function(q, scale) ifelse(q <= 0, 1, 2 * stats::pnorm(q, 0, scale, lower.tail = FALSE))
prior_phn = function(q, scale) {
  z = pmax(q, 0) / scale
  value = stats::pgamma(z^2 / 2, shape = .5)
  # Preserve tiny positive probabilities even when squaring z would underflow.
  ifelse(z < sqrt(.Machine$double.xmin), 2 * stats::dnorm(0) * z, value)
}
prior_qhn = function(p, scale) {
  value = stats::qnorm((1 - p) / 2, lower.tail = FALSE)
  small = which(p > 0 & p < sqrt(.Machine$double.eps))
  if (length(small)) {
    squared = stats::qgamma(p[small], shape = .5)
    value[small] = ifelse(squared == 0, p[small] / (2 * stats::dnorm(0)), sqrt(2) * sqrt(squared))
  }
  scale * value
}
prior_rhn = function(n, scale) abs(stats::rnorm(n, 0, scale))
prior_hn_mean = function(prior) prior$scale * sqrt(2 / pi)
prior_hn_sd = function(prior) prior$scale * sqrt(1 - 2 / pi)
prior_zero_mode = function(prior) 0
prior_undefined_moment = function(prior) NA_real_

prior_dht = function(x, scale, df) ifelse(x < 0, 0, 2 * stats::dt(x / scale, df) / scale)
prior_sht = function(q, scale, df) ifelse(q <= 0, 1, 2 * stats::pt(q / scale, df, lower.tail = FALSE))
prior_pht = function(q, scale, df) {
  z = pmax(q, 0) / scale
  log_ratio = 2 * log(z) - log(df)
  value = stats::pbeta(stats::plogis(log_ratio), .5, df / 2)
  ifelse(log_ratio < log(.Machine$double.xmin), 2 * stats::dt(0, df) * z, value)
}
prior_qht = function(p, scale, df) {
  value = stats::qt((1 - p) / 2, df, lower.tail = FALSE)
  small = which(p > 0 & p < sqrt(.Machine$double.eps))
  if (length(small)) {
    squared = stats::qbeta(p[small], .5, df / 2)
    value[small] = ifelse(squared == 0, p[small] / (2 * stats::dt(0, df)), sqrt(df) * sqrt(squared / (1 - squared)))
  }
  scale * value
}
prior_rht = function(n, scale, df) abs(scale * stats::rt(n, df))
prior_ht_mean = function(prior) {
  if (prior$df <= 1) return(NA_real_)
  # lbeta avoids cancellation between two large log-gamma values at large df.
  prior$scale * exp(log(2) + .5 * log(prior$df) - log(prior$df - 1) - lbeta(prior$df / 2, .5))
}
prior_ht_sd = function(prior) {
  if (prior$df <= 2) return(NA_real_)
  unit = prior_ht_mean(list(scale = 1, df = prior$df))
  prior$scale * sqrt(prior$df / (prior$df - 2) - unit^2)
}

prior_dhc = function(x, scale) prior_dht(x, scale, 1)
prior_shc = function(q, scale) prior_sht(q, scale, 1)
prior_phc = function(q, scale) prior_pht(q, scale, 1)
prior_qhc = function(p, scale) prior_qht(p, scale, 1)
prior_rhc = function(n, scale) prior_rht(n, scale, 1)
