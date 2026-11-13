test_that("prior checks reproduce exact registry draws and centered predictions", {
  p = prior_example()
  set.seed(91)
  a = check_prior(p, n = 20)
  set.seed(91)
  expect_equal(check_prior(p, n = 20), a)
  expect_s3_class(a, "flexbayes_prior_predictive")
  expect_equal(dim(a$parameter_draws), c(20L, 4L))
  expect_equal(dim(a$expected), c(20L, 6L))
  expect_equal(dim(a$simulated), c(20L, 6L))
  set.seed(91)
  for (name in names(p$selected)) {
    expect_equal(a$parameter_draws[, name], prior_distribution_operation(p$selected[[name]], "random", 20))
  }
  slopes = c("x", "gb")
  centered = sweep(p$model$X[, slopes], 2, p$model$reference, "-")
  expected = a$parameter_draws[, slopes] %*% t(centered) + a$parameter_draws[, "(Intercept)"]
  expect_equal(unname(a$expected), unname(expected))
  expect_equal(a$coefficient_draws[, "(Intercept)"], a$parameter_draws[, "(Intercept)"] - as.vector(a$parameter_draws[, slopes] %*% p$model$reference))
  noise = matrix(rnorm(120), 20)
  expect_equal(unname(a$simulated), unname(a$expected + noise * a$parameter_draws[, "sigma"]))
  expect_true(all(a$parameter_draws[, "sigma"] > 0))
  densities = plot(p)$data
  for (name in names(p$selected)) {
    rows = densities[densities$parameter == name, ]
    expect_equal(rows$density, prior_distribution_operation(p$selected[[name]], "density", rows$value))
  }
})

test_that("newdata preserves training centering, factor levels and outcome independence", {
  p = prior_example()
  new = data.frame(x = c(3.5, 10), g = c("a", "b"))
  set.seed(2)
  a = check_prior(p, n = 5, newdata = new)
  expect_equal(dim(a$expected), c(5L, 2L))
  expect_identical(colnames(a$X), colnames(p$model$X))
  expect_equal(a$metadata$reference, p$model$reference)
  new$y = c(1e10, -1e10)
  set.seed(2)
  expect_equal(check_prior(p, n = 5, newdata = new), a)
  expect_error(check_prior(p, newdata = data.frame(x = 1)), "missing predictors")
  expect_error(check_prior(p, newdata = data.frame(x = 1, g = "new")), "new level")
  expect_error(check_prior(p, newdata = data.frame(x = NA_real_, g = "a")))
  expect_error(check_prior(p, newdata = data.frame(x = Inf, g = "a")), "finite")
  single = identify_priors(y ~ x, p$model$data)
  draws = check_prior(single, 10, data.frame(x = mean(p$model$data$x)))
  expect_equal(as.vector(draws$expected), draws$parameter_draws[, "(Intercept)"])
  no_intercept = identify_priors(y ~ 0 + x, p$model$data)
  draws = check_prior(no_intercept, 2)
  expect_equal(draws$coefficient_draws, draws$parameter_draws[, "x", drop = FALSE])
  intercept = identify_priors(y ~ 1, p$model$data)
  draws = check_prior(intercept, 2)
  expect_equal(as.vector(draws$expected[, 1]), draws$parameter_draws[, "(Intercept)"])
})

test_that("t sampling uses selected scale and plots use the same mathematics", {
  p = prior_example()
  p = identify_priors(y ~ x, p$model$data)
  p = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  set.seed(8)
  result = check_prior(p, 30)
  set.seed(8)
  invisible(rnorm(30, p$selected$`(Intercept)`$mean, p$selected$`(Intercept)`$sd))
  expect_equal(result$parameter_draws[, "x"], rt(30, 7))
  densities = plot(p)$data
  rows = densities[densities$parameter == "x", ]
  expect_equal(rows$density, dt(rows$value, 7))
  expect_s3_class(plot(result) + ggplot2::theme_bw(), "ggplot")
  expect_equal(length(unique(plot(result)$data$response)), 2)
  expect_equal(unique(plot(result, type = "expected")$data$response), "Expected response")
  expect_equal(unique(plot(result, type = "simulated")$data$response), "Individual observation")
  expect_s3_class(plot(check_prior(p, 1, data.frame(x = 0))), "ggplot")
})

test_that("invalid priors and arguments fail without substitution", {
  p = prior_example()
  for (n in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) expect_error(check_prior(p, n), "positive integer")
  expect_error(check_prior(list()), "flexbayes_priors")
  q = update(p, parameter = "x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_error(check_prior(q), "contradictory")
  q = update(p, parameter = "x", belief(value < 0, probability = 0))
  expect_error(check_prior(q), "no_supported_fit")
  q = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  expect_error(check_prior(q), "backend_incompatible")
  q = p
  q$selected$x = NULL
  expect_error(check_prior(q), "missing")
  q$selected = p$selected
  q$selected$x$sd = Inf
  expect_error(check_prior(q), "Cannot sample prior for x")
  q$selected$x = list(distribution = "improper")
  expect_error(check_prior(q), "Unsupported distribution")
  q$selected$x = p$selected$x
  q$backend$compatible = FALSE
  expect_error(check_prior(q), "incompatible")
})
