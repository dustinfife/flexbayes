ratio_diagnostic_example = function(family = binomial()) {
  prior_example()
  identify_priors(y ~ x, data.frame(y = c(0, 1, 0, 1, 0, 1), x = (1:6) / 2), family)
}

test_that("ratio diagnostics disclose broad partial priors without changing them", {
  for (family in list(binomial(), poisson())) {
    p = ratio_diagnostic_example(family)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    p = (p + update_priors(parameter = "x", belief(mode = 2, scale = scale)))
    original = p
    compact = capture.output(print(p))
    text = gsub("[[:space:]]+", " ", paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n"))
    s = p$selected$x$sd
    mu = p$selected$x$mean
    expect_match(text, paste0("median = ", prior_ratio_format_log(mu)), fixed = TRUE)
    expect_match(text, prior_ratio_format_log(mu + qnorm(.025) * s), fixed = TRUE)
    expect_match(text, prior_ratio_format_log(mu + qnorm(.975) * s), fixed = TRUE)
    expect_match(text, paste0("Ratio mean = ", prior_ratio_format_log(mu + s^2 / 2)), fixed = TRUE)
    expect_match(text, "ratio mode = 2", fixed = TRUE)
    expect_match(text, "User-supplied constraint: ratio mode", fixed = TRUE)
    expect_match(text, "backend-default link SD retained", fixed = TRUE)
    expect_identical(p, original)
    expect_identical(capture.output(print(p)), compact)
    expect_false(any(grepl("Selected induced ratio", compact)))
    ordinary = capture.output(print(summary(p)))
    expect_false(any(grepl("Selected induced ratio", ordinary)))
    q = (ratio_diagnostic_example(family) + update_priors(parameter = "x", belief(sd = 2, scale = scale)))
    text = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
    expect_match(text, "backend-default link location retained")
  }
})

test_that("normal and t ratio summaries report valid quantiles and moment existence", {
  p = ratio_diagnostic_example()
  p = (p + update_priors(parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < exp(1), probability = pt(1, 7), scale = "odds_ratio")))
  expect_equal(p$selected$x$distribution, "student_t")
  text = gsub("[[:space:]]+", " ", paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n"))
  expect_match(text, "median = 1", fixed = TRUE)
  expect_match(text, prior_ratio_format_log(qt(.025, 7)), fixed = TRUE)
  expect_match(text, prior_ratio_format_log(qt(.975, 7)), fixed = TRUE)
  expect_match(text, "no finite mean")
  expect_match(text, "no attained positive global mode")
  expect_match(text, "no backend-default parameter was retained")
  p = (ratio_diagnostic_example() + update_priors(parameter = "x", belief(mean = 0), belief(sd = .5), belief(value < exp(.5), probability = pnorm(1), scale = "odds_ratio")))
  text = paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n")
  expect_match(text, paste0("Ratio mean = ", prior_ratio_format_log(.5^2 / 2)), fixed = TRUE)
  expect_match(text, paste0("ratio mode = ", prior_ratio_format_log(0 - .5^2)), fixed = TRUE)
})

test_that("extreme ratio summaries distinguish numerical range from undefined moments", {
  expect_equal(prior_ratio_format_log(log(245702)), "245,700")
  expect_match(prior_ratio_format_log(log(1e10)), "e\\+10")
  expect_match(prior_ratio_format_log(1000), "exceeds numeric range")
  expect_match(prior_ratio_format_log(-1000), "below numeric range")
  p = ratio_diagnostic_example()
  # Satisfiable probability constraints can accompany extreme link priors.
  p = (p + update_priors(parameter = "x", belief(mean = 1000), belief(sd = 1), belief(value < 1, probability = 0, scale = "odds_ratio")))
  # Use the diagnostic renderer directly on a selected valid mathematical prior
  # to exercise formatting independently of finite-tail fitting policy.
  p$selected$x = list(distribution = "normal", mean = 1000, sd = 1)
  text = paste(capture.output(prior_ratio_diagnostics(p, "x")), collapse = "\n")
  expect_match(text, "exceeds numeric range")
  expect_false(grepl("= Inf|= NaN", text))
})

test_that("quantile-reference diagnostics distinguish fixed and jointly resolved parameters", {
  p = ratio_diagnostic_example()
  p = (p + update_priors(parameter = "x", belief(value < 2, probability = .5, scale = "odds_ratio")))
  text = gsub("[[:space:]]+", " ", paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n"))
  expect_match(text, "User beliefs constrain link location; backend-default link SD retained")
  p = ratio_diagnostic_example()
  p = (p + update_priors(parameter = "x", belief(value < 2, probability = .8, scale = "odds_ratio")))
  text = gsub("[[:space:]]+", " ", paste(capture.output(print(summary(p, diagnostics = TRUE))), collapse = "\n"))
  expect_match(text, "neither link parameter was fixed")
})
