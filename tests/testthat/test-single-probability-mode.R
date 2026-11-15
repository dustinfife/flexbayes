single_mode_example = function() {
  prior_example()
  identify_priors(y ~ x, data.frame(y = c(0, 1, 0, 1, 0, 1), x = 1:6), binomial())
}

test_that("feasible single modes preserve actual backend spread", {
  for (target in c(.001, .999)) {
    p = single_mode_example()
    reference = p$defaults$`(Intercept)`
    expect_true(any(grepl("\\(Intercept\\).*Backend default", capture.output(print(p)))))
    q = update(p, parameter = "(Intercept)", belief(mode = target, scale = "probability"))
    expect_equal(q$selected$`(Intercept)`$sd, reference$sd)
    expect_equal(q$selected$`(Intercept)`$mean, qlogis(target) - reference$sd^2 * (2 * target - 1))
    expect_equal(prior_logit_normal_mode(q$selected$`(Intercept)`$mean, reference$sd), target, tolerance = 1e-8)
    expect_true(q$candidates$`(Intercept)`$reference_used)
    text = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
    expect_match(text, "Elicited*", fixed = TRUE)
    expect_match(gsub("[[:space:]]+", " ", text), "actual backend-default link SD = 2.5 preserved unchanged", fixed = TRUE)
    expect_equal(q$defaults, p$defaults)
  }
})

test_that("infeasible or competing global modes remain unresolved", {
  for (target in c(.5, .7, .95, .05)) {
    p = single_mode_example()
    q = update(p, parameter = "(Intercept)", belief(mode = target, scale = "probability"))
    expect_null(q$selected$`(Intercept)`)
    expect_identical(q$candidates$`(Intercept)`$status, "no_supported_fit")
    expect_match(q$candidates$`(Intercept)`$message, "Add an SD or informative quantile")
    expect_false(q$candidates$`(Intercept)`$alternatives$normal$reference_used)
    output = capture.output(print(q))
    expect_true(any(grepl("Unresolved", output)))
    expect_false(any(grepl("Elicited", output)))
  }
  expect_error(prior_logit_normal_mode(0, 2.5), "multiple global modes")
})

test_that("boundary modes and small-spread unimodality are checked", {
  p = single_mode_example()
  for (target in c(0, 1)) expect_error(update(p, parameter = "(Intercept)", belief(mode = target, scale = "probability")), "strictly inside")
  r = belief(mode = .7, scale = "probability")
  reference = list(distribution = "normal", mean = 0, sd = .8)
  fit = prior_fit_single_probability_mode(r, reference)
  expect_true(fit$valid)
  expect_equal(fit$distribution$sd, .8)
  transition = prior_fit_single_probability_mode(belief(mode = .5, scale = "probability"), list(distribution = "normal", mean = 0, sd = sqrt(2) * (1 - 1e-6)))
  expect_true(transition$valid)
  expect_equal(transition$achieved, .5)
  # No finite normal prior has its actual mode at the boundary; saturation
  # during numerical evaluation must not be accepted as a boundary mode.
  expect_equal(prior_logit_normal_mode(1000, 2.5), 1)
  extreme = prior_fit_single_probability_mode(belief(mode = 1e-300, scale = "probability"), p$defaults$`(Intercept)`)
  expect_true(is.logical(extreme$valid))
  saturated = prior_fit_single_probability_mode(belief(mode = .999, scale = "probability"), list(distribution = "normal", mean = 1000, sd = 100))
  expect_false(saturated$valid)
  expect_match(saturated$message, "no reference spread was modified")
  if (extreme$valid) expect_gt(extreme$achieved, 0)
})
