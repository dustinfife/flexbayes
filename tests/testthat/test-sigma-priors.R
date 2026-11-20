sigma_example = function() {
  skip_if_not_installed("rstanarm")
  if (as.character(packageVersion("rstanarm")) != "2.32.2") skip("Verified backend version required")
  identify_priors(sale_price ~ sqft + age + bedrooms, homes_2024)
}

test_that("a single sigma mean or SD preserves the default exponential family", {
  p = sigma_example()
  original = p
  for (type in c("mean", "sd")) {
    record = do.call(belief, stats::setNames(list(30000), type))
    q = p + update_priors("sigma", record)
    expect_identical(q$selected$sigma, list(distribution = "exponential", rate = 1 / 30000))
    expect_identical(q$defaults, p$defaults)
    expect_identical(q$selected[names(p$selected) != "sigma"], p$selected[names(p$selected) != "sigma"])
    expect_true(q$candidates$sigma$valid)
    expect_false(q$candidates$sigma$reference_used)
    expect_identical(q$beliefs$sigma[[1]]$type, type)
  }
  expect_identical(p, original)
  zero = p + update_priors("sigma", belief(mode = 0))
  expect_identical(zero$selected$sigma, p$defaults$sigma)
  expect_true(zero$candidates$sigma$reference_used)
  expect_match(paste(capture.output(print(zero)), collapse = "\n"), "Elicited*", fixed = TRUE)
})

test_that("positive-half distribution operations match independent formulas", {
  hn = list(distribution = "half_normal", scale = 2)
  ht = list(distribution = "half_student_t", scale = 2, df = 5)
  hc = list(distribution = "half_cauchy", scale = 2)
  expect_equal(prior_distribution_operation(hn, "mean"), 2 * sqrt(2 / pi))
  expect_equal(prior_distribution_operation(hn, "sd"), 2 * sqrt(1 - 2 / pi))
  m = 2 * sqrt(5 / pi) * gamma(2) / gamma(2.5)
  expect_equal(prior_distribution_operation(ht, "mean"), m)
  expect_equal(prior_distribution_operation(ht, "sd"), sqrt(4 * 5 / 3 - m^2))
  for (prior in list(hn, ht, hc)) {
    expect_equal(prior_distribution_operation(prior, "support"), c(0, Inf))
    expect_identical(prior_distribution_operation(prior, "mode"), 0)
    expect_equal(prior_distribution_operation(prior, "cdf", c(-1, 0)), c(0, 0))
    expect_equal(prior_distribution_operation(prior, "density", -1), 0)
    probability = c(.01, .5, .95, .999)
    quantiles = prior_distribution_operation(prior, "quantile", probability)
    expect_equal(prior_distribution_operation(prior, "cdf", quantiles), probability, tolerance = 1e-12)
    expect_equal(prior_distribution_operation(prior, "survival", quantiles), 1 - probability, tolerance = 1e-12)
    set.seed(2)
    expect_true(all(prior_distribution_operation(prior, "random", 100) > 0))
    arguments = c(list(f = prior_distribution_registry()[[prior$distribution]]$density, lower = 0, upper = Inf),
      prior[setdiff(names(prior), "distribution")])
    expect_equal(do.call(stats::integrate, arguments)$value, 1)
  }
  expect_true(is.na(prior_distribution_operation(hc, "mean")))
  expect_true(is.na(prior_distribution_operation(hc, "sd")))
  for (df in c(.5, 1, 1.5, 2)) {
    ht$df = df
    expect_true(is.na(prior_distribution_operation(ht, "sd")))
    expect_identical(is.na(prior_distribution_operation(ht, "mean")), df <= 1)
  }
})

test_that("sigma median and tail beliefs use the exponential rate parameterization", {
  p = sigma_example()
  median = p + update_priors("sigma", belief(value < 30000, probability = .5))
  expect_equal(median$selected$sigma$rate, log(2) / 30000)
  expect_equal(prior_distribution_operation(median$selected$sigma, "quantile", .5), 30000)
  for (operator in c("<", "<=", ">", ">=")) {
    record = belief(value < 60000, probability = .95)
    record$operator = operator
    q = p + update_priors("sigma", record)
    expected_rate = if (operator %in% c("<", "<=")) -log(.05) / 60000 else -log(.95) / 60000
    expect_equal(q$selected$sigma$rate, expected_rate)
    expect_lte(abs(prior_distribution_target(record, q$selected$sigma) - .95), record$tolerance)
  }
})

test_that("half-t candidates fit both df and scale without changing coefficient df policy", {
  p = sigma_example()
  q = p + update_priors("sigma", belief(mean = 30000), belief(value < 60000, probability = .95))
  expect_identical(q$selected$sigma$distribution, "half_student_t")
  expect_gt(q$selected$sigma$df, 1)
  expect_lt(q$selected$sigma$df, 2)
  expect_true(is.na(prior_distribution_operation(q$selected$sigma, "sd")))
  expect_lte(abs(prior_distribution_operation(q$selected$sigma, "mean") - 30000), 1e-6)
  expect_lte(abs(prior_distribution_operation(q$selected$sigma, "cdf", 60000) - .95), 1e-6)
  expect_identical(prior_student_t_policy()$df, 7)
  expect_identical(q$selected, (p + update_priors("sigma", belief(mean = 30000), belief(value < 60000, probability = .95)))$selected)
  expect_match(paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n"), "no finite moment", fixed = TRUE)
  # Recover several known half-t distributions, including an SD constraint.
  for (df in c(.7, 1.5, 3, 12)) {
    first = 100 * qt(.75, df)
    second = 100 * qt(.975, df)
    records = list(belief(value < first, probability = .5), belief(value < second, probability = .95))
    if (df > 2) records = list(belief(mean = prior_ht_mean(list(scale = 100, df = df))),
      belief(sd = prior_ht_sd(list(scale = 100, df = df))))
    result = prior_fit_scale_candidates(records, p$defaults$sigma)
    expect_true(result$valid)
    expect_identical(result$distribution$distribution, "half_student_t")
    expect_equal(result$distribution$df, df, tolerance = 1e-5)
    expect_equal(result$distribution$scale, 100, tolerance = 1e-5)
  }
})

test_that("sigma candidate ranks prefer defaults and simpler verified families", {
  p = sigma_example()
  exp = p + update_priors("sigma", belief(mean = 100), belief(sd = 100))
  expect_identical(exp$selected$sigma$distribution, "exponential")
  hn = p + update_priors("sigma", belief(mean = 100 * sqrt(2 / pi)), belief(sd = 100 * sqrt(1 - 2 / pi)))
  expect_identical(hn$selected$sigma$distribution, "half_normal")
  hc = p + update_priors("sigma", belief(value < 100, probability = .5),
    belief(value < 100 * qt(.975, 1), probability = .95))
  expect_identical(hc$selected$sigma$distribution, "half_cauchy")
  expect_equal(hc$candidates$sigma$rank, 1)
  expect_true(is.na(prior_distribution_operation(hc$selected$sigma, "mean")))
})

test_that("sigma failures retain beliefs and distinguish backend limits from contradictions", {
  p = sigma_example()
  mode = p + update_priors("sigma", belief(mode = 30000))
  expect_identical(mode$candidates$sigma$status, "backend_incompatible")
  expect_null(mode$selected$sigma)
  expect_equal(mode$beliefs$sigma[[1]]$target, 30000)
  expect_match(mode$candidates$sigma$message, "not reinterpreted", fixed = TRUE)
  expect_match(paste(capture.output(print(mode)), collapse = "\n"), "Incompatible")
  expect_error(check_prior(mode), "backend_incompatible")
  expect_error(as_rstanarm(mode), "Incompatible parameters: sigma")
  conflict = p + update_priors("sigma", belief(value < 0, probability = .5))
  expect_identical(conflict$candidates$sigma$status, "contradictory")
  domain = p + update_priors("sigma", belief(value <= 0, probability = 0))
  expect_identical(domain$selected$sigma, p$defaults$sigma)
  poor = p + update_priors("sigma", belief(mean = 30000), belief(sd = 1000))
  expect_null(poor$selected$sigma)
  expect_identical(poor$candidates$sigma$status, "no_supported_fit")
  expect_length(poor$beliefs$sigma, 2)
  expect_match(poor$candidates$sigma$message, "does not establish mathematical incompatibility", fixed = TRUE)
  expect_error(update(p, parameter = "sigma", belief(mean = .5, scale = "standardized")), "original")
  expect_identical(update(mode, parameter = "sigma", remove = 1)$selected$sigma, p$defaults$sigma)
})

test_that("selected sigma drives simulation, plotting and faithful prior_aux export", {
  p = sigma_example()
  variants = list(
    p + update_priors("sigma", belief(mean = 30000)),
    p + update_priors("sigma", belief(mean = 30000), belief(value < 60000, probability = .95)),
    p + update_priors("sigma", belief(mean = 100 * sqrt(2 / pi)), belief(sd = 100 * sqrt(1 - 2 / pi))),
    p + update_priors("sigma", belief(value < 100, probability = .5), belief(value < 100 * qt(.975, 1), probability = .95)))
  for (q in variants) {
    prior = q$selected$sigma
    native = as_rstanarm(q)$prior_aux
    expect_false(native$autoscale)
    expect_equal(native$scale, if (prior$distribution == "exponential") 1 / prior$rate else prior$scale)
    expect_equal(native$dist, switch(prior$distribution, exponential = "exponential", half_normal = "normal", "t"))
    if (prior$distribution == "half_student_t") expect_equal(native$df, prior$df)
    if (prior$distribution != "exponential") expect_equal(native$location, 0)
    # Capture backend preprocessing before sampling or Stan initialization.
    captured = new.env()
    testthat::local_mocked_bindings(set_sampling_args = function(..., data) {
      captured$data = data
      stop("captured sigma preprocessing")
    }, .package = "rstanarm")
    expect_error(do.call(rstanarm::stan_glm, as_rstanarm(q)), "captured sigma preprocessing")
    expect_equal(as.numeric(captured$data$prior_scale_for_aux), native$scale)
    expect_equal(as.numeric(captured$data$prior_dist_for_aux), switch(prior$distribution, exponential = 3, half_normal = 1, 2))
    if (prior$distribution == "half_student_t") expect_equal(as.numeric(captured$data$prior_df_for_aux), prior$df)
    set.seed(45)
    check = check_prior(q, n = 10000, newdata = q$model$data[1:2, ])
    set.seed(45)
    expected_draws = prior_sample_parameters(q, 10000)
    noise = matrix(rnorm(length(check$expected)), nrow = nrow(check$expected))
    expect_identical(check$parameter_draws, expected_draws)
    expect_equal(check$simulated, check$expected + noise * expected_draws[, "sigma"])
    expect_equal(as.numeric(quantile(check$parameter_draws[, "sigma"], c(.5, .95))),
      prior_distribution_operation(prior, "quantile", c(.5, .95)), tolerance = .12)
    plot = plot(q, parameter = "sigma") + ggplot2::theme_bw()
    expect_s3_class(plot, "ggplot")
    selected = plot$data[grepl("selected", plot$data$Candidate, fixed = TRUE), ]
    expect_equal(selected$density, prior_distribution_operation(prior, "density", selected$value))
  }
})

test_that("outcome reconciliation never fits or sacrifices the sigma prior", {
  p = joint_example() + update_priors("sigma", belief(mean = 4))
  before = p$selected$sigma
  q = p + update_outcome(list(x = 6, z = 1), belief(mean = 12))
  expect_identical(q$selected$sigma, before)
  expect_identical(q$beliefs$sigma, p$beliefs$sigma)
  expect_true(q$joint$status %in% c("compatible", "reconciled"))
  expect_no_error(check_prior(q, 10))
  revised = update(q, parameter = "sigma", remove = 1, belief(mean = 2))
  expect_identical(revised$selected$sigma, list(distribution = "exponential", rate = .5))
  expect_equal(revised$selected[names(revised$selected) != "sigma"], q$selected[names(q$selected) != "sigma"])
  expect_no_error(as_rstanarm(revised))
})

test_that("sigma fitting uses acceptance bands and preserves measurement units", {
  p = sigma_example()
  approximate = p + update_priors("sigma", belief(mean = 100, tolerance = 10), belief(sd = 110))
  expect_identical(approximate$selected$sigma$distribution, "exponential")
  expect_lte(abs(1 / approximate$selected$sigma$rate - 100), 10)
  expect_lte(abs(1 / approximate$selected$sigma$rate - 110), 1e-6)
  expect_equal(approximate$beliefs$sigma[[1]]$target, 100)
  expect_equal(approximate$beliefs$sigma[[1]]$tolerance, 10)
  for (probability in c(0, 1)) {
    record = belief(value < 10, probability = probability)
    candidate = p + update_priors("sigma", record)
    expect_true(candidate$candidates$sigma$valid)
    expect_lte(abs(prior_distribution_target(record, candidate$selected$sigma) - probability), 1e-6)
    expect_identical(candidate$beliefs$sigma[[1]]$target, probability)
  }
  q = p + update_priors("sigma", belief(mean = 30000), belief(value < 60000, probability = .95))
  d = p$model$data
  d$sale_price = d$sale_price / 1000
  units = identify_priors(sale_price ~ sqft + age + bedrooms, d) +
    update_priors("sigma", belief(mean = 30, tolerance = 1e-9), belief(value < 60, probability = .95))
  expect_equal(units$selected$sigma$df, q$selected$sigma$df, tolerance = 1e-8)
  expect_equal(units$selected$sigma$scale * 1000, q$selected$sigma$scale, tolerance = 1e-7)
})

test_that("sigma numerical failures do not become mathematical contradictions", {
  p = sigma_example()
  local_mocked_bindings(prior_fit_half_t = function(...) stop("forced numerical failure"))
  q = p + update_priors("sigma", belief(mean = 30000), belief(value < 60000, probability = .95))
  expect_identical(q$candidates$sigma$status, "no_supported_fit")
  expect_length(q$beliefs$sigma, 2)
  expect_null(q$selected$sigma)
  expect_match(q$candidates$sigma$message, "forced numerical failure", fixed = TRUE)
  expect_error(check_prior(q), "no_supported_fit")
})

test_that("fit uses the elicited prior_aux without requiring routine Stan sampling", {
  p = sigma_example() + update_priors("sigma", belief(mean = 30000))
  captured = new.env()
  local_mocked_bindings(prior_rstanarm_fit = function(arguments) {
    captured$arguments = arguments
    "fit sentinel"
  })
  expect_identical(fit(p, chains = 1, iter = 20), "fit sentinel")
  expect_equal(captured$arguments$prior_aux, as_rstanarm(p)$prior_aux)
  expect_equal(captured$arguments$prior_aux$scale, 30000)
  expect_equal(captured$arguments$chains, 1)
})

test_that("half-t large-df moments are stable and sigma does not change expected responses", {
  prior = list(distribution = "half_student_t", scale = 100, df = 1e12)
  expect_equal(prior_distribution_operation(prior, "mean"), 100 * sqrt(2 / pi), tolerance = 1e-9)
  expect_equal(prior_distribution_operation(prior, "sd"), 100 * sqrt(1 - 2 / pi), tolerance = 1e-9)
  expect_error(prior_validate_distribution(list(distribution = "half_student_t", scale = 1, df = 0)), "positive")
  p = sigma_example()
  q = p + update_priors("sigma", belief(mean = 30000))
  set.seed(13)
  before = check_prior(p, 100)
  set.seed(13)
  after = check_prior(q, 100)
  expect_identical(before$expected, after$expected)
  expect_false(identical(before$simulated, after$simulated))
  expect_equal(mean(after$parameter_draws[, "sigma"]), 30000, tolerance = .2)
})

test_that("joint summaries retain sigma beliefs and non-exponential sigma selections", {
  p = joint_example() + update_priors("sigma", belief(mean = 4), belief(value < 8, probability = .95))
  before = p$selected$sigma
  q = p + update_outcome(list(x = 6, z = 1), belief(mean = 12))
  expect_identical(q$selected$sigma, before)
  expect_identical(q$beliefs$sigma, p$beliefs$sigma)
  text = paste(capture.output(print(summary(q))), collapse = "\n")
  expect_match(text, "Residual-scale prior: Half-t", fixed = TRUE)
  expect_match(text, "1. mean = 4", fixed = TRUE)
  expect_match(text, "no finite moment", fixed = TRUE)
  expect_no_error(as_rstanarm(q))
  expect_no_error(check_prior(q, 100))
})

test_that("positive-half tails preserve tiny probabilities without cancellation", {
  for (family in c("half_normal", "half_student_t", "half_cauchy")) {
    prior = prior_scale_distribution(family, 1, if (family == "half_student_t") 5 else NULL)
    probabilities = c(1e-200, 1e-20, 1e-10)
    quantiles = prior_distribution_operation(prior, "quantile", probabilities)
    expect_true(all(quantiles > 0))
    expect_equal(prior_distribution_operation(prior, "cdf", quantiles) / probabilities,
      rep(1, 3), tolerance = 1e-10)
    expect_gt(prior_distribution_operation(prior, "cdf", 1e-20), 0)
  }
  p = sigma_example()
  p$selected$sigma = list(distribution = "exponential", rate = 1e-309)
  expect_error(plot(p, parameter = "sigma"), "numerical plotting range")
})
