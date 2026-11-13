glm_prior_example = function(family) {
  prior_example()
  d = data.frame(y = c(0, 1, 0, 1, 1, 0), x = 1:6, g = factor(rep(c("a", "b"), 3)))
  if (family$family == "poisson") d$y = c(0, 2, 1, 4, 0, 3)
  identify_priors(y ~ x + g, d, family)
}

test_that("GLM effective defaults match actual backend preprocessing", {
  captured = new.env()
  testthat::local_mocked_bindings(set_sampling_args = function(..., data) {
    captured$data = data
    stop("captured GLM")
  }, .package = "rstanarm")
  for (family in list(binomial(), poisson())) {
    p = glm_prior_example(family)
    expect_false("sigma" %in% names(p$selected))
    expect_equal(p$defaults$`(Intercept)`, list(distribution = "normal", mean = 0, sd = 2.5))
    expect_equal(p$defaults$x$sd, 2.5 / sd(1:6))
    expect_error(rstanarm::stan_glm(y ~ x + g, data = p$model$data, family = family), "captured GLM")
    expect_equal(as.numeric(captured$data$prior_scale), c(p$defaults$x$sd, p$defaults$gb$sd))
    expect_equal(as.numeric(captured$data$prior_mean_for_intercept), 0)
    expect_equal(as.numeric(captured$data$prior_scale_for_intercept), 2.5)
    expect_equal(as.numeric(captured$data$xbar), unname(p$model$reference))
    a = as_rstanarm(p)
    expect_null(a$prior_aux)
    expect_false("prior_aux" %in% names(a))
    expect_error(do.call(rstanarm::stan_glm, c(a, list(prior_PD = TRUE))), "captured GLM")
    expect_false(a$prior$autoscale)
    expect_error(do.call(rstanarm::stan_glm, a), "captured GLM")
    expect_equal(as.numeric(captured$data$prior_scale), c(p$defaults$x$sd, p$defaults$gb$sd))
    zero = p$model$data
    zero$y = 0
    z = identify_priors(y ~ x + g, zero, family)
    expect_error(do.call(rstanarm::stan_glm, as_rstanarm(z)), "captured GLM")
  }
})

test_that("GLM prediction inverse links and sampling are exact and reproducible", {
  for (family in list(binomial(), poisson())) {
    p = glm_prior_example(family)
    # Moderate priors avoid the deliberately guarded Poisson extreme tail.
    p = update(p, parameter = "(Intercept)", belief(mean = 0), belief(sd = .2))
    p = update(p, parameter = "x", belief(mean = .1), belief(sd = .1))
    p = update(p, parameter = "gb", belief(mean = 0), belief(sd = .1))
    set.seed(51)
    a = check_prior(p, 10)
    set.seed(51)
    expect_equal(check_prior(p, 10), a)
    expect_equal(dim(a$expected), c(10L, 6L))
    expect_equal(ncol(a$parameter_draws), 3L)
    centered = sweep(p$model$X[, c("x", "gb")], 2, p$model$reference, "-")
    eta = a$parameter_draws[, c("x", "gb")] %*% t(centered) + a$parameter_draws[, "(Intercept)"]
    expect_equal(unname(a$linear_predictor), unname(eta))
    expected = if (family$family == "binomial") plogis(eta) else exp(eta)
    expect_equal(unname(a$expected), unname(expected))
    set.seed(51)
    for (name in names(p$selected)) invisible(prior_distribution_operation(p$selected[[name]], "random", 10))
    values = if (family$family == "binomial") rbinom(60, 1, expected) else rpois(60, expected)
    expect_equal(as.vector(a$simulated), values)
    expect_true(all(a$simulated >= 0 & a$simulated == floor(a$simulated)))
    if (family$family == "binomial") expect_true(all(a$expected >= 0 & a$expected <= 1 & a$simulated <= 1))
    b = check_prior(p, 4, data.frame(x = 3.5, g = "a"))
    expect_equal(dim(b$expected), c(4L, 1L))
    expect_s3_class(plot(a) + ggplot2::theme_bw(), "ggplot")
  }
})

test_that("GLM beliefs retain link scale and compatibility rules", {
  for (family in list(binomial(), poisson())) {
    p = glm_prior_example(family)
    q = update(p, parameter = "x", belief(mode = .5), belief(sd = .2))
    expect_equal(q$selected$x$mean, .5)
    expect_equal(unname(as_rstanarm(q)$prior$location[1]), .5)
    expect_error(update(p, parameter = "x", belief(mean = 1, scale = "standardized")), "Gaussian-style full standardization")
    expect_error(plot(p, parameter = "x", scale = "standardized"), "Gaussian-style")
    t = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
    expect_identical(t$candidates$x$status, "backend_incompatible")
    expect_error(as_rstanarm(t), "Incompatible parameters")
    t = update(t, parameter = "gb", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
    expect_equal(t$selected$x$distribution, "student_t")
    expect_equal(as_rstanarm(t)$prior$dist, "t")
    bad = update(p, parameter = "x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
    expect_error(check_prior(bad), "contradictory")
    expect_error(as_rstanarm(bad), "contradictory")
  }
})

test_that("GLM responses, links and numerical limits are validated", {
  p = glm_prior_example(binomial())
  d = p$model$data
  for (family in list(binomial("probit"), poisson("identity"), Gamma())) expect_error(identify_priors(y ~ x, d, family), "Supported families")
  d$y = rep(.5, 6)
  expect_error(identify_priors(y ~ x, d, binomial()), "binary 0/1")
  expect_error(identify_priors(cbind(y, y) ~ x, d, binomial()), "response|Response")
  d$y = factor(rep(c("a", "b"), 3))
  expect_error(identify_priors(y ~ x, d, binomial()), "response|Response")
  for (y in list(c(-1, 0:4), c(.5, 0:4), rep(Inf, 6))) {
    d$y = y
    expect_error(identify_priors(y ~ x, d, poisson()))
  }
  d$y = rep(0, 6)
  expect_s3_class(identify_priors(y ~ x, d, poisson()), "flexbayes_priors")
  expect_s3_class(identify_priors(y ~ x, d, binomial()), "flexbayes_priors")
  d$y = rep(2, 6)
  expect_error(identify_priors(y ~ x, d, poisson()), "constant responses")
  eta = matrix(c(-1e300, 1e300), 1)
  a = prior_simulate_response(eta, matrix(numeric(), 1, 0), binomial())
  expect_equal(as.vector(a$expected), c(0, 1))
  expect_error(prior_simulate_response(matrix(1000), matrix(numeric(), 1, 0), poisson()), "Poisson predicted rates")
})

test_that("fit forwards GLM family and returns native backend result", {
  captured = new.env()
  native = structure(list(), class = c("stanreg", "glm", "lm"))
  testthat::local_mocked_bindings(prior_rstanarm_fit = function(arguments) { captured$args = arguments; native })
  for (family in list(binomial(), poisson())) {
    p = glm_prior_example(family)
    expect_identical(fit(p, chains = 2), native)
    expect_equal(captured$args$family$family, family$family)
    expect_equal(captured$args$family$link, family$link)
    expect_null(captured$args$prior_aux)
  }
})

test_that("optional GLM fits return stanreg models with explicit priors", {
  skip_if(Sys.getenv("FLEXBAYES_RUN_STAN") != "true", "Enable optional Stan sampling")
  for (family in list(binomial(), poisson())) {
    p = glm_prior_example(family)
    m = fit(p, chains = 2, iter = 2000, warmup = 1000, seed = 506, refresh = 0)
    expect_s3_class(m, "stanreg")
    expect_equal(m$family$family, family$family)
    expect_false("sigma" %in% colnames(as.matrix(m)))
    expect_equal(unname(m$prior.info$prior$scale), c(p$selected$x$sd, p$selected$gb$sd))
  }
})
