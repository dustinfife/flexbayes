outcome_example = function(family = binomial()) {
  prior_example()
  identify_priors(y ~ x, data.frame(y = c(0, 1, 0, 1, 0, 1), x = 1:6), family)
}

test_that("transformed normal moments and modes are mathematically correct", {
  prior = list(distribution = "normal", mean = .5, sd = .4)
  for (type in c("mean", "sd", "mode")) {
    r = do.call(belief, c(stats::setNames(list(.5), type), list(scale = "response")))
    expected = switch(type, mean = exp(.5 + .4^2 / 2), sd = exp(.5 + .4^2 / 2) * sqrt(expm1(.4^2)), mode = exp(.5 - .4^2))
    expect_equal(prior_distribution_target(r, prior), expected)
  }
  set.seed(61)
  mc = plogis(rnorm(300000, .5, .4))
  expect_lt(abs(prior_distribution_target(belief(mean = .6, scale = "probability"), prior) - mean(mc)), .001)
  expect_lt(abs(prior_distribution_target(belief(sd = .1, scale = "probability"), prior) - sd(mc)), .001)
  mode = prior_distribution_target(belief(mode = .6, scale = "probability"), prior)
  expect_equal(prior_logit_mode_equation(qlogis(mode), .5, .4^2), 0, tolerance = 1e-9)
  expect_gt(abs(mode - plogis(.5)), .001)
  expect_error(prior_logit_normal_mode(0, 2), "multiple global modes")
})

test_that("outcome quantiles mix with link beliefs and boundaries", {
  p = outcome_example()
  p = update(p, parameter = "(Intercept)", belief(mean = 0), belief(value < plogis(1), probability = pnorm(1), scale = "probability"))
  expect_equal(p$selected$`(Intercept)`$sd, 1, tolerance = 1e-8)
  q = update(p, parameter = "(Intercept)", belief(value > plogis(2), probability = .95, scale = "probability"))
  expect_identical(q$candidates$`(Intercept)`$status, "contradictory")
  q = update(outcome_example(), parameter = "(Intercept)", belief(value < 0, probability = .1, scale = "probability"))
  expect_identical(q$candidates$`(Intercept)`$status, "contradictory")
  q = update(outcome_example(), parameter = "(Intercept)", belief(value <= 1, probability = 1, scale = "probability"))
  expect_equal(q$selected$`(Intercept)`, q$defaults$`(Intercept)`)
  q = update(outcome_example(), parameter = "(Intercept)", belief(value < .5, probability = 0, scale = "probability"))
  expect_null(q$selected$`(Intercept)`)
  p = outcome_example(poisson())
  p = update(p, parameter = "(Intercept)", belief(mean = log(3)), belief(value <= 6, probability = .95, scale = "response"))
  expect_equal(p$selected$`(Intercept)`$sd, log(2) / qnorm(.95), tolerance = 1e-8)
})

test_that("nonlinear outcome modes and moments fit and export link priors", {
  for (family in list(binomial(), poisson())) {
    p = outcome_example(family)
    scale = if (family$family == "binomial") "probability" else "response"
    target = if (scale == "probability") .7 else 3
    threshold = if (scale == "probability") .9 else 6
    p = update(p, parameter = "(Intercept)", do.call(belief, list(mode = target, scale = scale)),
      belief(value < threshold, probability = .95, scale = scale))
    expect_equal(p$candidates$`(Intercept)`$status, "selected")
    expect_equal(p$candidates$`(Intercept)`$achieved, c(target, .95), tolerance = 1e-6)
    expect_equal(as_rstanarm(p)$prior_intercept$dist, "normal")
    expect_equal(p$selected$`(Intercept)`$distribution, "normal")
    plot = plot(p, parameter = "(Intercept)", scale = scale)
    values = plot$data$value
    link = prior_outcome_link(values, scale)
    jac = if (scale == "probability") values * (1 - values) else values
    prior = p$selected$`(Intercept)`
    expect_equal(plot$data$density, dnorm(link, prior$mean, prior$sd) / jac)
    expect_s3_class(plot + ggplot2::theme_bw(), "ggplot")
    text = paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n")
    expect_match(text, "centering reference")
    expect_match(text, "Induced candidate families")
    expect_match(text, "Outcome moments")
    set.seed(83)
    check = check_prior(p, n = 5, newdata = data.frame(x = mean(p$model$data$x)))
    alpha = check$parameter_draws[, "(Intercept)"]
    expect_equal(as.vector(check$expected), if (scale == "probability") plogis(alpha) else exp(alpha))
    arguments = as_rstanarm(p)
    expect_equal(unname(arguments$prior_intercept$location), prior$mean)
    expect_equal(unname(arguments$prior_intercept$scale), prior$sd)
  }
})

test_that("nonlinear transformed means and SDs fit without naive transforms", {
  p = outcome_example(poisson())
  q = update(p, parameter = "(Intercept)", belief(mean = 3, scale = "response"), belief(sd = 1, scale = "response"))
  expected_sd = sqrt(log1p(1 / 9))
  expect_equal(q$selected$`(Intercept)`$mean, log(3) - expected_sd^2 / 2, tolerance = 1e-6)
  expect_equal(q$selected$`(Intercept)`$sd, expected_sd, tolerance = 1e-6)
  expect_false(q$candidates$`(Intercept)`$alternatives$student_t$valid)
  expect_match(q$candidates$`(Intercept)`$alternatives$student_t$message, "no finite positive moments")
  p = outcome_example()
  q = update(p, parameter = "(Intercept)", belief(mean = .3), belief(sd = .1, scale = "probability"))
  expect_equal(q$selected$`(Intercept)`$mean, .3, tolerance = 1e-6)
  expect_equal(prior_distribution_target(q$beliefs$`(Intercept)`[[2]], q$selected$`(Intercept)`), .1, tolerance = 1e-6)
  q = update(p, parameter = "(Intercept)", belief(mean = .7, scale = "probability"))
  expect_true(q$candidates$`(Intercept)`$reference_used)
  expect_equal(q$candidates$`(Intercept)`$achieved, .7, tolerance = 1e-6)
  expect_true(is.finite(prior_distribution_target(belief(sd = .1, scale = "probability"), list(distribution = "student_t", location = 0, scale = 1, df = 7))))
})

test_that("outcome domains, unsupported modes and removal are explicit", {
  p = outcome_example()
  expect_error(update(p, parameter = "x", belief(mean = .7, scale = "probability")), "intercept only")
  expect_error(update(p, parameter = "(Intercept)", belief(mean = 1, scale = "probability")), "strictly inside")
  expect_error(update(p, parameter = "(Intercept)", belief(sd = .5, scale = "probability")), "below 0.5")
  expect_error(update(p, parameter = "(Intercept)", belief(value < -1, probability = .5, scale = "probability")), "in ")
  expect_error(update(p, parameter = "(Intercept)", belief(mean = 3, scale = "response")), "family")
  q = update(p, parameter = "(Intercept)", belief(mode = .7, scale = "probability"))
  expect_equal(q$candidates$`(Intercept)`$achieved, .7, tolerance = 1e-6)
  expect_match(q$candidates$`(Intercept)`$alternatives$student_t$message, "global mode")
  q = update(q, parameter = "(Intercept)", remove = 1)
  expect_equal(q$selected$`(Intercept)`, p$defaults$`(Intercept)`)
  p = outcome_example(poisson())
  expect_error(update(p, parameter = "(Intercept)", belief(mean = 0, scale = "response")), "positive")
  expect_error(update(p, parameter = "(Intercept)", belief(value < -1, probability = .5, scale = "response")), "positive")
  q = update(outcome_example(), parameter = "(Intercept)", belief(mean = .7, scale = "probability", tolerance = 1e-10))
  expect_null(q$selected$`(Intercept)`)
  expect_match(q$candidates$`(Intercept)`$message, "numerical evaluation")
})

test_that("redundant outcome beliefs preserve reference anchoring and diagnostics", {
  p = outcome_example()
  p = update(p, parameter = "(Intercept)", belief(mean = .7, scale = "probability"), belief(mean = .7, scale = "probability"))
  expect_true(p$candidates$`(Intercept)`$reference_used)
  expect_length(p$candidates$`(Intercept)`$achieved, 2)
  expect_equal(p$selected$`(Intercept)`$sd, p$defaults$`(Intercept)`$sd)
  expect_silent(capture.output(print(summary(p, diagnostics = TRUE))))
})

test_that("general nonlinear rank deficiency remains explicit", {
  p = outcome_example()
  p = update(p, parameter = "(Intercept)", belief(mean = 0), belief(mean = .5, scale = "probability"))
  expect_null(p$selected$`(Intercept)`)
  expect_match(p$candidates$`(Intercept)`$message, "Underidentified nonlinear")
})

test_that("fitting outcome beliefs forwards only selected link-scale priors", {
  captured = new.env()
  testthat::local_mocked_bindings(prior_rstanarm_fit = function(arguments) { captured$args = arguments; structure(list(), class = "stanreg") })
  for (family in list(binomial(), poisson())) {
    p = outcome_example(family)
    original = p$model$data
    scale = if (family$family == "binomial") "probability" else "response"
    target = if (scale == "probability") .7 else 3
    p = update(p, parameter = "(Intercept)", do.call(belief, list(mean = target, scale = scale)))
    expect_s3_class(fit(p, chains = 2), "stanreg")
    expect_equal(captured$args$data, original)
    expect_equal(unname(captured$args$prior_intercept$location), p$selected$`(Intercept)`$mean)
    expect_equal(unname(captured$args$prior_intercept$scale), p$selected$`(Intercept)`$sd)
    expect_equal(captured$args$family$family, family$family)
  }
})
