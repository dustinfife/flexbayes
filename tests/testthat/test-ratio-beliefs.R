ratio_example = function(family = binomial()) {
  prior_example()
  identify_priors(y ~ x, data.frame(y = c(0, 1, 0, 1, 0, 1), x = 1:6), family)
}

test_that("ratio quantiles exactly match link constraints", {
  for (family in list(binomial(), poisson())) {
    p = ratio_example(family)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    for (ratio in c(.5, 1, 2)) {
      a = update(p, parameter = "x", belief(value < ratio, probability = .8, scale = scale))
      b = update(p, parameter = "x", belief(value < log(ratio), probability = .8))
      expect_equal(a$selected$x, b$selected$x)
      a = update(p, parameter = "x", belief(value > ratio, probability = .8, scale = scale))
      b = update(p, parameter = "x", belief(value > log(ratio), probability = .8))
      expect_equal(a$selected$x, b$selected$x)
      expect_identical(a$beliefs$x[[1]]$scale, scale)
      expect_equal(a$beliefs$x[[1]]$threshold, ratio)
    }
  }
})

test_that("single ratio moments retain actual backend reference parameters", {
  for (family in list(binomial(), poisson())) {
    p = ratio_example(family)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    s = p$defaults$x$sd
    for (target in c(.5, 1, 2)) {
      mode = update(p, parameter = "x", belief(mode = target, scale = scale))
      expect_equal(mode$selected$x$mean, log(target) + s^2)
      expect_equal(mode$selected$x$sd, s)
      mean = update(p, parameter = "x", belief(mean = target, scale = scale))
      expect_equal(mean$selected$x$mean, log(target) - s^2 / 2)
      expect_equal(mean$selected$x$sd, s)
      sd = update(p, parameter = "x", belief(sd = target, scale = scale))
      expect_equal(sd$selected$x$mean, p$defaults$x$mean)
      expect_equal(prior_distribution_target(sd$beliefs$x[[1]], sd$selected$x), target)
      expect_true(mode$candidates$x$reference_used)
      expect_false(mode$candidates$x$alternatives$student_t$valid)
      expect_match(mode$candidates$x$alternatives$student_t$message, "no finite positive moments")
    }
  }
})

test_that("ratio multi-beliefs, density, export and native fitting boundary agree", {
  captured = new.env()
  testthat::local_mocked_bindings(prior_rstanarm_fit = function(arguments) { captured$args = arguments; structure(list(), class = "stanreg") })
  for (family in list(binomial(), poisson())) {
    p = ratio_example(family)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    q = update(p, parameter = "x", belief(mode = 2, scale = scale), belief(value < 4, probability = .95, scale = scale))
    expect_equal(q$candidates$x$status, "selected")
    expect_equal(q$candidates$x$achieved, c(2, .95), tolerance = 1e-6)
    plot = plot(q, parameter = "x", scale = scale)
    expect_s3_class(plot + ggplot2::theme_bw(), "ggplot")
    expect_equal(plot$data$density, dnorm(log(plot$data$value), q$selected$x$mean, q$selected$x$sd) / plot$data$value)
    expect_true(all(plot$data$value > 0))
    expect_equal(plot$scales$get_scales("x")$trans$name, "log-10")
    expect_s3_class(fit(q, iter = 100), "stanreg")
    expect_equal(captured$args$data, p$model$data)
    expect_equal(unname(captured$args$prior$location), q$selected$x$mean)
    expect_equal(unname(as_rstanarm(q)$prior$scale), q$selected$x$sd)
    check = check_prior(q, n = 3, newdata = data.frame(x = mean(p$model$data$x)))
    expect_equal(as.vector(check$expected), if (family$family == "binomial") plogis(check$parameter_draws[, "(Intercept)"]) else exp(check$parameter_draws[, "(Intercept)"]))
    text = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
    expect_match(text, scale)
    expect_match(gsub("[[:space:]]+", " ", text), "one-unit predictor increase")
    expect_match(text, "lognormal")
    expect_match(text, "Target.*Achieved")
    mixed = update(p, parameter = "x", belief(mean = log(2)), belief(value < 4, probability = .95, scale = scale))
    expect_equal(mixed$selected$x$sd, log(2) / qnorm(.95), tolerance = 1e-8)
    moments = update(p, parameter = "x", belief(mean = 2, scale = scale), belief(sd = 1, scale = scale))
    expect_equal(moments$selected$x$sd, sqrt(log1p(1 / 4)), tolerance = 1e-6)
    expect_equal(moments$selected$x$mean, log(2) - log1p(1 / 4) / 2, tolerance = 1e-6)
  }
})

test_that("invalid ratio scales, targets and terms are rejected", {
  p = ratio_example()
  b = belief(mode = 2, scale = "odds_ratio")
  expect_error(update(p, parameter = "(Intercept)", b), "numeric slope")
  expect_error(update(p, parameter = "x", belief(mode = 2, scale = "rate_ratio")), "family")
  expect_error(update(ratio_example(poisson()), parameter = "x", b), "family")
  expect_error(update(prior_example(), parameter = "x", b), "family")
  for (value in c(0, -1)) {
    expect_error(update(p, parameter = "x", belief(mode = value, scale = "odds_ratio")), "positive finite")
    expect_error(update(p, parameter = "x", belief(value < value, probability = .5, scale = "odds_ratio")), "positive finite")
  }
  for (value in c(NA, Inf, NaN)) expect_error(belief(mean = value, scale = "odds_ratio"))
  d = p$model$data
  d$g = factor(rep(c("a", "b"), 3))
  for (formula in list(y ~ x + I(x^2), y ~ x + log(x), y ~ log(x), y ~ poly(x, 2), y ~ x * g, y ~ I(x^2), y ~ g)) {
    q = identify_priors(formula, d, binomial())
    name = setdiff(names(q$selected), "(Intercept)")[1]
    expect_error(update(q, parameter = name, b), "untransformed numeric slope")
  }
})

test_that("ratio contradictions and t-only compatibility remain distinct", {
  p = ratio_example()
  q = update(p, parameter = "x", belief(value < 1, probability = .95, scale = "odds_ratio"), belief(value > 2, probability = .95, scale = "odds_ratio"))
  expect_identical(q$candidates$x$status, "contradictory")
  d = p$model$data
  d$z = c(1, 0, 3, 0, 5, 2)
  p = identify_priors(y ~ x + z, d, binomial())
  q = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < exp(1), probability = pt(1, 7), scale = "odds_ratio"))
  expect_identical(q$candidates$x$status, "backend_incompatible")
  expect_equal(q$candidates$x$alternatives$student_t$distribution$df, 7)
  expect_error(as_rstanarm(q), "backend_incompatible.*allow_fallback")
})

test_that("ratio plot range and retained-reference diagnostics are explicit", {
  p = ratio_example()
  q = update(p, parameter = "x", belief(mode = 2, scale = "odds_ratio"))
  text = gsub("[[:space:]]+", " ", paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n"))
  expect_match(text, "Elicited*", fixed = TRUE)
  expect_match(text, "Backend-default link SD")
  expect_match(text, "retained unchanged; location elicited", fixed = TRUE)
  q = update(p, parameter = "x", belief(sd = 2, scale = "odds_ratio"))
  expect_match(q$candidates$x$numerical, "Backend-default link mean = 0 retained unchanged", fixed = TRUE)
  q = update(p, parameter = "x", belief(mean = 1000), belief(sd = 1))
  expect_error(plot(q, parameter = "x", scale = "odds_ratio"), "no density points were silently removed")
})
