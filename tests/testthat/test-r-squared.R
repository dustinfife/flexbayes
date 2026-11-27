r2_example = function(package = "rstanarm") {
  skip_if_not_installed(package)
  identify_priors(y ~ x, data.frame(x = seq(-1, 1, length.out = 12), y = seq(-2, 2, length.out = 12)), package = package)
}
r2_exp_probability = function(h, coefficient_sd, predictor_sd, rate) {
  a = sqrt((1 - h) / h) * coefficient_sd * predictor_sd * rate
  1 - 2 * exp(a^2 / 2 + pnorm(-a, log.p = TRUE))
}
r2_fixed = function(package = "rstanarm") {
  p = r2_example(package)
  p + update_priors("x", belief(mean = 2), belief(sd = 1e-10, tolerance = 1e-14))
}

test_that("R-squared targets and probability-only semantics validate", {
  p = r2_example()
  for (value in c(-.1, 1.1)) expect_error(p + update_priors("sigma",
    belief(value < value, probability = .5, scale = "r_squared")))
  for (b in list(belief(mean = .7, scale = "r_squared"), belief(mode = .7, scale = "r_squared"), belief(sd = .1, scale = "r_squared"))) {
    expect_error(p + update_priors("sigma", b), "probability/quantile")
  }
  expect_error(p + update_priors("x", belief(value < .7, probability = .5, scale = "r_squared")), "sigma in a Gaussian")
  p = identify_priors(y ~ x, data.frame(x = 1:6, y = rep(0:1, 3)), family = binomial())
  expect_error(p + update_priors("x", belief(value < .7, probability = .5, scale = "r_squared")), "Gaussian")
})

test_that("sigma is fitted using uncertain coefficients and exact analytic probability", {
  p = r2_example()
  before = p
  set.seed(123)
  seed = .Random.seed
  q = p + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  expect_identical(.Random.seed, seed)
  expect_identical(p, before)
  expect_identical(q$selected$x, p$selected$x)
  expect_identical(q$selected$`(Intercept)`, p$selected$`(Intercept)`)
  expect_identical(q$r_squared$status, "verified")
  expect_identical(q$selected$sigma$distribution, "exponential")
  actual = r2_exp_probability(.6, p$selected$x$sd, sd(p$model$data$x), q$selected$sigma$rate)
  expect_lte(abs(actual - .9), 1e-6)
  info = q$r_squared$verification
  expect_lte(abs(actual - info$achieved), info$error + 1e-12)
  expect_lte(abs(info$achieved - .9) + info$error, 1e-6)
  expect_true(info$precise)
  expect_identical(q, p + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared")))
  median = p + update_priors("sigma", belief(value < .7, probability = .5, scale = "r_squared"))
  expect_lte(abs(1 - r2_exp_probability(.7, p$selected$x$sd, sd(p$model$data$x), median$selected$sigma$rate) - .5), 1e-6)
  a = as_rstanarm(q)
  expect_equal(a$prior_aux$scale, 1 / q$selected$sigma$rate)
  set.seed(5)
  check = check_prior(q, n = 20000)
  v = check$parameter_draws[, "x"]^2 * var(p$model$data$x)
  r2 = v / (v + check$parameter_draws[, "sigma"]^2)
  expect_equal(mean(r2 > .6), .9, tolerance = .015)
  expect_equal(check$expected, check$coefficient_draws %*% t(check$X))
})

test_that("covariance uses all estimation model-matrix columns and their correlations", {
  p = identify_priors(mpg ~ wt * hp + I(wt^2) + factor(am), mtcars)
  expect_equal(p$model$predictor_covariance, cov(p$model$X[, -1, drop = FALSE]))
  p = identify_priors(y ~ x + z, data.frame(x = 1:12, z = (1:12) + rep(c(-1, 1), 6), y = seq_len(12)))
  context = prior_r_squared_context(p)
  grid = prior_r_squared_grid(context, 64)
  expected = sum(diag(context$covariance) * vapply(p$selected[c("x", "z")], function(s) s$sd^2, numeric(1)))
  expect_equal(sum(grid$v * grid$w), expected, tolerance = 1e-7)
  diagonal = context
  diagonal$B = diag(diag(context$B))
  g = prior_r_squared_grid(diagonal, 64)
  prior = list(distribution = "exponential", rate = .3)
  b = belief(value > .5, probability = .8, scale = "r_squared")
  expect_gt(abs(prior_r_squared_probability(b, prior, grid) - prior_r_squared_probability(b, prior, g)), .001)
  saved = p$model$predictor_covariance
  q = p + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  r = q + update_priors("x", belief(mode = 1))
  expect_identical(r$model$predictor_covariance, saved)
  expect_identical(r$r_squared$context$coefficients, r$selected[c("x", "z")])
  stale = q
  stale$model$data$x[1] = -100
  expect_error(check_prior(stale), "changed")
  stale = q
  stale$model$predictor_covariance[1, 1] = -1
  expect_error(check_prior(stale), "covariance changed")
})

test_that("multiple R-squared beliefs preserve sigma family ranking", {
  p = r2_example()
  s = p$selected$x$sd * sd(p$model$data$x)
  a = r2_exp_probability(.5, p$selected$x$sd, sd(p$model$data$x), .8)
  b = 1 - r2_exp_probability(.85, p$selected$x$sd, sd(p$model$data$x), .8)
  q = p + update_priors("sigma", belief(value > .5, probability = a, scale = "r_squared"),
    belief(value < .85, probability = b, scale = "r_squared"))
  expect_identical(q$selected$sigma$distribution, "exponential")
  expect_equal(q$selected$sigma$rate, .8, tolerance = 1e-6)
  # Independent Normal/half-Normal ratio has an exact arctangent probability.
  a = 2 / pi * atan(s * sqrt((1 - .5) / .5) / 1.2)
  b = 1 - 2 / pi * atan(s * sqrt((1 - .85) / .85) / 1.2)
  q = p + update_priors("sigma", belief(value > .5, probability = a, scale = "r_squared"),
    belief(value < .85, probability = b, scale = "r_squared"))
  expect_identical(q$selected$sigma$distribution, "half_normal")
  expect_lte(abs(2 / pi * atan(s / q$selected$sigma$scale) - a), 1e-6)
  expect_lte(abs(1 - 2 / pi * atan(s * sqrt(.15 / .85) / q$selected$sigma$scale) - b), 1e-6)
  expect_true(q$r_squared$verification$valid)
})

test_that("effectively fixed coefficients use a verified uncertainty envelope", {
  p = r2_fixed()
  v = 4 * var(p$model$data$x)
  probability = pexp(sqrt(v * .4 / .6), rate = .7)
  q = p + update_priors("sigma", belief(value > .6, probability = probability, scale = "r_squared"))
  expect_equal(q$selected$sigma$rate, .7, tolerance = 1e-7)
  context = prior_r_squared_context(q)
  grid = prior_r_squared_grid(context, 64)
  expect_length(grid$bounds, 2)
  expect_lte(grid$bounds[1], v)
  expect_gte(grid$bounds[2], v)
  expect_lte(abs(q$r_squared$verification$achieved - probability), q$r_squared$verification$error + 1e-9)
})

test_that("direct and derived sigma beliefs are evaluated jointly and can switch backend", {
  skip_if_not_installed("brms")
  p = r2_fixed()
  v = 4 * var(p$model$data$x)
  probability = pgamma(sqrt(v * .4 / .6), shape = 3, rate = 1)
  q = p + update_priors("sigma", belief(mode = 2),
    belief(value > .6, probability = probability, scale = "r_squared"))
  expect_identical(q$backend$package, "brms")
  expect_identical(q$selected$sigma$distribution, "gamma")
  expect_false(q$candidates$sigma$reference_used)
  expect_lte(abs(prior_distribution_operation(q$selected$sigma, "mode") - 2), 1e-6)
  expect_lte(abs(pgamma(sqrt(v * .4 / .6), q$selected$sigma$shape, q$selected$sigma$rate) - probability), 1e-6)
  expect_s3_class(as_brms(q)$prior, "brmsprior")
  r = p + update_priors("sigma", belief(mean = 1 / .7),
    belief(value > .6, probability = pexp(sqrt(v * .4 / .6), .7), scale = "r_squared"))
  expect_identical(r$selected$sigma$distribution, "exponential")
  expect_length(r$beliefs$sigma, 2)
})

test_that("coefficient changes refit sigma without changing the requested coefficients", {
  p = r2_example()
  q = p + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  requested = p + update_priors("x", belief(mode = 2))
  r = q + update_priors("x", belief(mode = 2))
  expect_identical(r$selected$x, requested$selected$x)
  expect_identical(r$r_squared$status, "verified")
  expect_false(identical(r$selected$sigma, q$selected$sigma))
  expect_identical(r$beliefs$sigma, q$beliefs$sigma)
  expect_identical(r$r_squared$context$coefficients, r$selected["x"])
  # Intercept changes have no effect on the R-squared distribution.
  i = q + update_priors("(Intercept)", belief(mode = 100))
  expect_identical(i$selected$sigma, q$selected$sigma)
  expect_identical(i$r_squared, q$r_squared)
  ordered = requested + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  expect_equal(r$selected, ordered$selected, tolerance = 1e-7)
  stale = q
  stale$selected$x$sd = 100
  expect_error(check_prior(stale), "stale")
  removed = q + update_priors("sigma", remove = 1)
  expect_null(removed$r_squared)
  expect_identical(removed$selected$sigma, p$defaults$sigma)
})

test_that("joint outcomes and differences precede derived sigma elicitation", {
  p = r2_example()
  difference = update_priors(difference = list(from = list(x = -1), to = list(x = 1)), belief(mean = 4))
  outcome = update_priors(outcome = list(x = 0), belief(mode = 3))
  rs = update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  q = p + difference + outcome + rs
  r = p + rs + outcome + difference
  expect_identical(q$r_squared$status, "verified")
  expect_true(all(q$joint$outcome$satisfied))
  expect_identical(q$r_squared$context$coefficients, q$selected["x"])
  expect_equal(q$selected, r$selected, tolerance = 1e-6)
  expect_silent(prior_joint_validate(q))
  expect_silent(prior_r_squared_validate(q))
  expect_s3_class(check_prior(q, n = 3), "flexbayes_prior_predictive")
})

test_that("failed dependency refits retain a coherent collection with explicit safeguards", {
  q = r2_example() + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  saved = q$selected
  local_mocked_bindings(prior_fit_scale_candidates = function(...) list(valid = FALSE,
    status = "no_supported_fit", alternatives = list(), message = "Forced numerical verification failure."))
  expect_warning({ r = q + update_priors("x", belief(mode = 20)) }, "numerical verification")
  expect_identical(r$selected, saved)
  expect_identical(r$r_squared$status, "unresolved")
  expect_length(r$beliefs$x, 1)
  expect_identical(r$beliefs$sigma, q$beliefs$sigma)
  expect_error(check_prior(r), "Unsuccessful")
  expect_warning(check_prior(r, n = 3, allow_fallback = TRUE), "fallback")
  expect_warning(as_rstanarm(r, allow_fallback = TRUE), "fallback")
})

test_that("mathematical impossibility is distinct from numerical precision limits", {
  p = r2_example()
  expect_warning({ q = p + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"),
    belief(value < .5, probability = .9, scale = "r_squared")) }, "contradict|probabilit")
  expect_identical(q$candidates$sigma$status, "contradictory")
  for (b in list(belief(value < 0, probability = .1, scale = "r_squared"),
    belief(value > 1, probability = .1, scale = "r_squared"))) {
    expect_warning({ q = p + update_priors("sigma", b) }, "mathematically incompatible")
    expect_identical(q$candidates$sigma$status, "contradictory")
  }
  expect_warning({ q = p + update_priors("sigma", belief(value > .6, probability = .9,
    scale = "r_squared", tolerance = 1e-12)) }, "verified|unresolved")
  expect_identical(q$r_squared$status, "unresolved")
  expect_identical(q$beliefs$sigma[[1]]$tolerance, 1e-12)
  expect_identical(q$selected, p$selected)
  zero = identify_priors(y ~ 1, data.frame(y = 1:8))
  q = zero + update_priors("sigma", belief(value <= 0, probability = 1, scale = "r_squared"))
  expect_identical(q$r_squared$status, "verified")
  expect_identical(q$selected$sigma, zero$selected$sigma)
  expect_warning({ q = zero + update_priors("sigma", belief(value > 0, probability = .1, scale = "r_squared")) }, "mathematically")
})

test_that("summaries and plots distinguish derived R-squared from marginal sigma", {
  q = r2_example() + update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  text = paste(capture.output(print(summary(q))), collapse = "\n")
  expect_match(text, "Beliefs about model R-squared")
  expect_match(text, "empirical predictor distribution")
  expect_match(text, "Requested.*Achieved")
  detailed = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
  expect_match(detailed, "Numerical_error")
  expect_match(detailed, "covariance")
  expect_match(detailed, "shortcut")
  set.seed(2)
  saved = .Random.seed
  plot = plot(q, parameter = "sigma", scale = "r_squared")
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot + ggplot2::theme_minimal()))
  expect_identical(.Random.seed, saved)
  expect_true(all(diff(plot$data$probability) >= -1e-10))
  expect_s3_class(plot(q, parameter = "sigma"), "ggplot")
  expect_error(plot(q, parameter = "x", scale = "r_squared"), "sigma")
})

test_that("Student's t coefficient uncertainty uses its actual density", {
  p = r2_example() + update_priors("x", belief(mean = 0), belief(sd = sqrt(7 / 5)),
    belief(value < 1, probability = pt(1, 7)))
  expect_identical(p$selected$x$distribution, "student_t")
  q = p + update_priors("sigma", belief(value > .6, probability = .8, scale = "r_squared"))
  sx = sd(p$model$data$x)
  rate = q$selected$sigma$rate
  exact = integrate(function(z) (1 - exp(-abs(z) * sx * sqrt(.4 / .6) * rate)) * dt(z, 7),
    -Inf, Inf, abs.tol = 1e-10, rel.tol = 1e-10)$value
  expect_lte(abs(exact - .8), 1e-6)
  expect_lte(abs(exact - q$r_squared$verification$achieved), q$r_squared$verification$error + 1e-9)
})

test_that("three-dimensional quadrature matches an independent chi distribution", {
  d = expand.grid(x = c(-1, 1), z = c(-1, 1), w = c(-1, 1))
  d$y = 1:8
  p = identify_priors(y ~ x + z + w, d)
  q = p + update_priors("sigma", belief(value > .6, probability = .8, scale = "r_squared"))
  expect_identical(q$r_squared$status, "verified")
  s = p$selected$x$sd * sd(d$x)
  rate = q$selected$sigma$rate
  exact = integrate(function(r) (1 - exp(-s * r * sqrt(.4 / .6) * rate)) *
    2 * r * dchisq(r^2, 3), 0, Inf, abs.tol = 1e-10, rel.tol = 1e-10)$value
  expect_lte(abs(exact - .8), 1e-6)
  expect_lte(abs(exact - q$r_squared$verification$achieved), q$r_squared$verification$error + 1e-9)
})

test_that("bounded higher-dimensional searches remain unresolved without dropping columns", {
  d = as.data.frame(expand.grid(a = c(-1, 1), b = c(-1, 1), c = c(-1, 1), d = c(-1, 1)))
  d$y = seq_len(nrow(d))
  p = identify_priors(y ~ ., d)
  expect_warning({ q = p + update_priors("sigma", belief(value > .6, probability = .8, scale = "r_squared")) }, "three varying|unresolved|verified")
  expect_identical(q$r_squared$status, "unresolved")
  expect_identical(q$selected, p$selected)
  expect_length(q$r_squared$context$coefficients, 4)
  boundary = p + update_priors("sigma", belief(value > 0, probability = 1, scale = "r_squared", tolerance = 1e-14))
  expect_identical(boundary$r_squared$status, "verified")
  expect_identical(boundary$selected$sigma, p$selected$sigma)
})

test_that("equivalent backends and unit changes preserve induced probabilities", {
  skip_if_not_installed("brms")
  p = r2_example()
  r = update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  q = p + r
  b = r2_example("brms") + r
  expect_equal(q$selected, b$selected)
  expect_equal(q$r_squared$verification$achieved, b$r_squared$verification$achieved)
  expect_s3_class(as_brms(b)$prior, "brmsprior")
  d = p$model$data
  d$x = d$x * 1000
  u = identify_priors(y ~ x, d) + r
  expect_equal(q$selected$sigma, u$selected$sigma, tolerance = 1e-7)
  restored = switch_package(q, "brms")
  expect_equal(restored$selected, q$selected)
})

test_that("mathematical outcome incompatibility survives a blocked sigma dependency", {
  q = r2_example() + update_priors(outcome = list(x = 0), belief(mean = 1)) +
    update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  expect_warning({ bad = q + update_priors(outcome = list(x = 0), belief(mean = 10)) }, "incompatible|verified")
  expect_identical(bad$joint$status, "incompatible")
  expect_warning(check_prior(bad, n = 2, allow_fallback = TRUE), "fallback")
})

test_that("half-t and lognormal sigma candidates fit actual induced quantiles", {
  p = r2_fixed()
  v = 4 * var(p$model$data$x)
  a = prior_pht(sqrt(v * .7 / .3), 1.2, 5)
  b = prior_pht(sqrt(v * .2 / .8), 1.2, 5)
  q = p + update_priors("sigma", belief(value > .3, probability = a, scale = "r_squared"),
    belief(value > .8, probability = b, scale = "r_squared"))
  expect_identical(q$selected$sigma$distribution, "half_student_t")
  expect_equal(q$selected$sigma$scale, 1.2, tolerance = 1e-4)
  expect_equal(q$selected$sigma$df, 5, tolerance = 1e-4)
  expect_lte(abs(prior_pht(sqrt(v * .7 / .3), q$selected$sigma$scale, q$selected$sigma$df) - a), 1e-6)
  skip_if_not_installed("brms")
  sdlog = sqrt(log(1.5) / 1.5)
  meanlog = log(2) + sdlog^2
  probability = plnorm(sqrt(v * .4 / .6), meanlog, sdlog)
  q = p + update_priors("sigma", belief(mode = 2), belief(mean = 3),
    belief(value > .6, probability = probability, scale = "r_squared"))
  expect_identical(q$backend$package, "brms")
  expect_identical(q$selected$sigma$distribution, "lognormal")
  expect_equal(q$selected$sigma$meanlog, meanlog, tolerance = 1e-6)
  expect_equal(q$selected$sigma$sdlog, sdlog, tolerance = 1e-6)
  expect_false(q$candidates$sigma$reference_used)
})

test_that("observed outcome variance is diagnostic rather than an extra constraint", {
  p = r2_example()
  record = belief(value > .6, probability = .9, scale = "r_squared")
  context = prior_r_squared_context(p)
  a = prior_fit_scale_candidates(list(record), p$defaults$sigma, context = context)
  context$observed_variance = 1e12
  b = prior_fit_scale_candidates(list(record), p$defaults$sigma, context = context)
  expect_identical(a$distribution, b$distribution)
  expect_identical(a$achieved, b$achieved)
})

test_that("mixed-scale contradiction numbers refer to current sigma positions", {
  p = r2_example()
  expect_warning({ q = p + update_priors("sigma", belief(mean = 1),
    belief(value > .6, probability = .9, scale = "r_squared"),
    belief(value < .5, probability = .9, scale = "r_squared")) }, "Beliefs 2 and 3")
  expect_identical(q$candidates$sigma$status, "contradictory")
})

test_that("the housing model retains coefficient priors while eliciting sigma", {
  p = identify_priors(sale_price ~ sqft + age + bedrooms, homes_2024)
  q = p + update_priors("sigma", belief(value < .7, probability = .5, scale = "r_squared"))
  expect_identical(q$r_squared$status, "verified")
  expect_identical(q$selected[1:4], p$selected[1:4])
  expect_lte(abs(q$r_squared$verification$achieved - .5) + q$r_squared$verification$error, 1e-6)
  expect_identical(q$model$predictor_covariance, p$model$predictor_covariance)
  expect_equal(as_rstanarm(q)$prior_aux$scale, 1 / q$selected$sigma$rate)
  set.seed(321)
  expect_s3_class(check_prior(q, n = 20), "flexbayes_prior_predictive")
})

test_that("overlapping probability bands preserve every entered answer", {
  q = r2_example() + update_priors("sigma",
    belief(value > .6, probability = .9, scale = "r_squared"),
    belief(value > .6, probability = .9000015, scale = "r_squared"))
  expect_identical(q$r_squared$status, "verified")
  expect_identical(vapply(q$beliefs$sigma, `[[`, numeric(1), "target"), c(.9, .9000015))
  expect_true(all(abs(q$r_squared$verification$achieved - c(.9, .9000015)) + q$r_squared$verification$error <= 1e-6))
})
