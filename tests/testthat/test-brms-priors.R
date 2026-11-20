brms_example = function(family = gaussian(), formula = y ~ x + z) {
  skip_if_not_installed("brms")
  d = data.frame(y = c(2, 4, 3, 7, 6, 10), x = 1:6, z = rep(0:1, 3), g = factor(rep(c("a", "b"), 3)))
  if (family$family == "binomial") d$y = c(0, 0, 1, 0, 1, 1)
  identify_priors(formula, d, family, package = "brms")
}

brms_t_beliefs = function(p, parameter) update(p, parameter = parameter,
  belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))

test_that("proper default policy is independent of fitting backend", {
  skip_if_not_installed("rstanarm")
  for (family in list(gaussian(), binomial(), poisson())) {
    b = brms_example(family)
    r = identify_priors(b$model$formula, b$model$data, family)
    expect_identical(b$defaults, r$defaults)
    expect_identical(b$default_origin, r$default_origin)
    expect_identical(b$model$reference, r$model$reference)
    expect_identical(b$selected, r$selected)
    set.seed(24)
    a = check_prior(b, 30, newdata = data.frame(x = 3.5, z = .5))
    set.seed(24)
    z = check_prior(r, 30, newdata = data.frame(x = 3.5, z = .5))
    expect_identical(a$parameter_draws, z$parameter_draws)
    expect_identical(a$expected, z$expected)
    expect_identical(a$simulated, z$simulated)
    text = paste(capture.output(print(summary(b))), collapse = "\n")
    expect_match(gsub("\\s+", " ", text), "native brms defaults are overridden", fixed = TRUE)
    expect_match(text, "flexbayes default", fixed = TRUE)
  }
})

test_that("brms mapping and centering match real design matrices and Stan code", {
  for (formula in list(y ~ x * g, y ~ x + z, y ~ 0 + x + g, y ~ I(x^2) + g)) {
    p = brms_example(formula = formula)
    args = as_brms(p)
    sd = do.call(brms::make_standata, args)
    code = do.call(brms::make_stancode, args)
    expect_equal(unname(sd$X), unname(p$model$X), ignore_attr = TRUE)
    slopes = setdiff(colnames(p$model$X), "(Intercept)")
    expect_setequal(args$prior$coef[args$prior$class == "b"], unname(p$backend$coefficients))
    expect_length(p$backend$coefficients, length(slopes))
    expect_true(all(nzchar(args$prior$prior)))
    if ("(Intercept)" %in% colnames(p$model$X)) {
      expect_equal(p$model$reference, colMeans(p$model$X[, slopes, drop = FALSE]))
      expect_match(code, "means_X[i - 1] = mean(X[, i])", fixed = TRUE)
      expect_match(code, "Intercept - dot_product(means_X, b)", fixed = TRUE)
      expect_equal(sum(args$prior$class == "Intercept"), 1)
    } else expect_false("Intercept" %in% args$prior$class)
  }
  p = brms_example(formula = y ~ x * g)
  contrasts(p$model$data$g) = contr.sum(2)
  p = identify_priors(y ~ x * g, p$model$data, package = "brms")
  args = as_brms(p)
  expect_equal(contrasts(args$data$g), contrasts(p$model$data$g))
  expect_equal(unname(do.call(brms::make_standata, args)$X), unname(p$model$X), ignore_attr = TRUE)
  altered = p
  altered$model$reference[1] = 10
  expect_error(as_brms(altered), "Centering reference")
  altered = p
  altered$model$data$x[1] = 100
  expect_error(as_brms(altered), "changed")
  altered = p
  altered$backend$coefficients[1] = "wrong"
  expect_error(as_brms(altered), "mapping changed")
})

test_that("brms supports coefficient-specific families and binary/count exports", {
  p = brms_t_beliefs(brms_example(), "x")
  expect_identical(p$selected$x$distribution, "student_t")
  expect_identical(p$selected$z, p$defaults$z)
  expect_true(p$backend$compatible)
  expect_length(p$fallbacks, 0)
  args = as_brms(p)
  expect_match(args$prior$prior[args$prior$coef == "x" & args$prior$class == "b"], "student_t(7,", fixed = TRUE)
  expect_match(args$prior$prior[args$prior$coef == "z" & args$prior$class == "b"], "normal(", fixed = TRUE)
  expect_s3_class(plot(p, parameter = "x"), "ggplot")
  expect_s3_class(check_prior(p, 20), "flexbayes_prior_predictive")
  for (family in list(binomial(), poisson())) {
    q = brms_example(family)
    args = as_brms(q)
    expect_equal(args$family$family, if (family$family == "binomial") "bernoulli" else "poisson")
    expect_false("sigma" %in% args$prior$class)
    expect_type(do.call(brms::make_stancode, args), "character")
  }
})

test_that("Gamma and lognormal registry operations use actual positive moments", {
  gamma = list(distribution = "gamma", shape = 4, rate = .5)
  ln = list(distribution = "lognormal", meanlog = 2, sdlog = .4)
  for (p in list(gamma, ln)) {
    expect_equal(prior_distribution_operation(p, "support"), c(0, Inf))
    q = prior_distribution_operation(p, "quantile", c(.1, .5, .95))
    expect_equal(prior_distribution_operation(p, "cdf", q), c(.1, .5, .95))
    expect_equal(prior_distribution_operation(p, "survival", q), c(.9, .5, .05))
    expect_true(all(prior_distribution_operation(p, "random", 10) > 0))
    expect_true(all(prior_distribution_operation(p, "density", q) > 0))
  }
  expect_equal(prior_distribution_operation(gamma, "mean"), 8)
  expect_equal(prior_distribution_operation(gamma, "sd"), 4)
  expect_equal(prior_distribution_operation(gamma, "mode"), 6)
  expect_equal(prior_distribution_operation(ln, "mean"), exp(2 + .4^2 / 2))
  expect_equal(prior_distribution_operation(ln, "sd"), exp(2 + .4^2 / 2) * sqrt(expm1(.4^2)))
  expect_equal(prior_distribution_operation(ln, "mode"), exp(2 - .4^2))
  expect_error(prior_validate_distribution(list(distribution = "gamma", shape = 0, rate = 1)), "positive")
  expect_error(prior_validate_distribution(list(distribution = "lognormal", meanlog = 0, sdlog = -1)), "positive")
})

test_that("positive-mode sigma beliefs retain a disclosed reference SD", {
  skip_if_not_installed("rstanarm")
  p = brms_example()
  q = p + update_priors("sigma", belief(mode = 3))
  expect_identical(q$selected$sigma$distribution, "gamma")
  expect_equal(prior_distribution_operation(q$selected$sigma, "mode"), 3, tolerance = 1e-6)
  expect_equal(prior_distribution_operation(q$selected$sigma, "sd"), 1 / p$defaults$sigma$rate, tolerance = 1e-6)
  expect_true(q$candidates$sigma$reference_used)
  expect_match(q$candidates$sigma$numerical, "default's sd", fixed = TRUE)
  expect_length(q$fallbacks, 0)
  expect_s3_class(plot(q, parameter = "sigma"), "ggplot")
  args = as_brms(q)
  expect_match(args$prior$prior[args$prior$class == "sigma"], "gamma(", fixed = TRUE)
  code = do.call(brms::make_stancode, args)
  expect_match(code, "gamma_lpdf(sigma", fixed = TRUE)
  expect_match(code, "real<lower=0> sigma", fixed = TRUE)
  expect_match(paste(capture.output(print(summary(q))), collapse = "\n"), "mode: 3", fixed = TRUE)
  r = identify_priors(p$model$formula, p$model$data) + update_priors("sigma", belief(mode = 3))
  expect_identical(r$candidates$sigma$status, "backend_incompatible")
})

test_that("multiple positive sigma constraints fit without dropping beliefs", {
  p = brms_example()
  for (family in c("gamma", "lognormal")) {
    target = if (family == "gamma") list(distribution = "gamma", shape = 4, rate = .3) else list(distribution = "lognormal", meanlog = 2, sdlog = .7)
    mode = prior_distribution_operation(target, "mode")
    mean = prior_distribution_operation(target, "mean")
    threshold = prior_distribution_operation(target, "quantile", .95)
    q = p + update_priors("sigma", belief(mode = mode), belief(mean = mean), belief(value < threshold, probability = .95))
    expect_identical(q$selected$sigma$distribution, family)
    expect_length(q$beliefs$sigma, 3)
    expect_false(q$candidates$sigma$reference_used)
    expect_equal(prior_distribution_operation(q$selected$sigma, "mode"), mode, tolerance = 1e-6)
    expect_equal(prior_distribution_operation(q$selected$sigma, "mean"), mean, tolerance = 1e-6)
    expect_equal(prior_distribution_operation(q$selected$sigma, "cdf", threshold), .95, tolerance = 1e-6)
    code = do.call(brms::make_stancode, as_brms(q))
    expect_match(code, paste0(family, "_lpdf(sigma"), fixed = TRUE)
    set.seed(8)
    expected = prior_distribution_operation(q$selected$sigma, "random", 500)
    set.seed(8)
    # Sample sigma directly with the same public registry operation used in simulation.
    expect_identical(expected, prior_distribution_operation(q$selected$sigma, "random", 500))
    check = check_prior(q, 1000)
    expect_equal(mean(check$parameter_draws[, "sigma"]), mean, tolerance = .15)
    density = plot(q, parameter = "sigma")$data
    expect_true(all(is.finite(density$density)))
    selected = density[grepl("(selected)", density$Candidate, fixed = TRUE), ]
    expect_equal(selected$density, prior_distribution_operation(q$selected$sigma, "density", selected$value))
    expect_identical(as_brms(q)$prior$prior[as_brms(q)$prior$class == "sigma"], prior_brms_distribution(q$selected$sigma))
  }
})

test_that("brms preserves fallback, recovery, and outcome reconciliation", {
  skip_if_not_installed("rstanarm")
  p = brms_example() + update_priors("sigma", belief(mean = 3))
  q = p + update_priors("sigma", belief(value < 0, probability = .95))
  expect_identical(q$selected, p$selected)
  expect_identical(q$candidates$sigma$status, "contradictory")
  expect_error(as_brms(q), "allow_fallback")
  expect_error(check_prior(q), "allow_fallback")
  expect_error(fit(q), "allow_fallback")
  expect_warning({args = as_brms(q, allow_fallback = TRUE)}, "fallback")
  expect_equal(args, as_brms(p))
  expect_warning(check_prior(q, 5, allow_fallback = TRUE), "fallback")
  expect_s3_class(plot(q), "ggplot")
  expect_null(update(q, parameter = "sigma", remove = 2)$fallbacks$sigma)
  b = brms_example() + update_outcome(list(x = 3, z = 0), belief(mean = 4), belief(value < 4 - qnorm(.95), probability = .05))
  r = identify_priors(b$model$formula, b$model$data) + update_outcome(list(x = 3, z = 0), belief(mean = 4), belief(value < 4 - qnorm(.95), probability = .05))
  expect_identical(b$selected, r$selected)
  expect_identical(b$selected$sigma, b$defaults$sigma)
  expect_type(do.call(brms::make_stancode, as_brms(b)), "character")
})

test_that("brms fit dispatch protects parameterization and consumes fallback approval", {
  p = brms_example()
  captured = new.env()
  local_mocked_bindings(prior_brms_fit = function(arguments) { captured$arguments = arguments; "sentinel" })
  expect_identical(fit(p, chains = 1, sample_prior = "only"), "sentinel")
  expect_identical(captured$arguments$sample_prior, "only")
  expect_error(fit(p, prior = "normal(0,1)"), "Cannot override")
  q = p + update_priors("x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_warning(expect_identical(fit(q, allow_fallback = TRUE, chains = 1), "sentinel"), "fallback")
  expect_false("allow_fallback" %in% names(captured$arguments))
})

test_that("brms native prior-only smoke fit uses exported proper priors", {
  skip_if(Sys.getenv("FLEXBAYES_TEST_STAN") != "true", "Set FLEXBAYES_TEST_STAN=true for compilation and sampling")
  p = brms_example(formula = y ~ x)
  fit = fit(p, chains = 1, iter = 200, warmup = 100, seed = 12, refresh = 0, sample_prior = "only")
  expect_s3_class(fit, "brmsfit")
  expect_true(all(is.finite(as.matrix(fit))))
})

test_that("brms preserves transformed belief scales and numeric sigma failures", {
  skip_if_not_installed("rstanarm")
  p = brms_example()
  r = identify_priors(p$model$formula, p$model$data)
  b = p + update_priors("x", belief(mean = 1, scale = "standardized"), belief(sd = .5, scale = "standardized"))
  q = r + update_priors("x", belief(mean = 1, scale = "standardized"), belief(sd = .5, scale = "standardized"))
  expect_identical(b$selected, q$selected)
  for (family in list(binomial(), poisson())) {
    p = brms_example(family)
    r = identify_priors(p$model$formula, p$model$data, family)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    constraints = list(belief(value < 2, probability = .5, scale = scale), belief(value < 4, probability = .95, scale = scale))
    b = do.call(update, c(list(p, parameter = "x"), constraints))
    q = do.call(update, c(list(r, parameter = "x"), constraints))
    expect_identical(b$selected, q$selected)
    expect_type(do.call(brms::make_stancode, as_brms(b)), "character")
    scale = if (family$family == "binomial") "probability" else "response"
    threshold = if (family$family == "binomial") .7 else 3
    b = b + update_priors("(Intercept)", belief(value < threshold, probability = .8, scale = scale))
    q = q + update_priors("(Intercept)", belief(value < threshold, probability = .8, scale = scale))
    expect_identical(b$selected, q$selected)
    expect_type(do.call(brms::make_stancode, as_brms(b)), "character")
  }
  p = brms_example() + update_priors("sigma", belief(mode = 3))
  q = p + update_priors("sigma", belief(mean = 1))
  expect_identical(q$selected, p$selected)
  expect_identical(q$candidates$sigma$status, "no_supported_fit")
  expect_match(q$candidates$sigma$message, "does not establish mathematical incompatibility", fixed = TRUE)
  expect_warning(as_brms(q, allow_fallback = TRUE), "fallback")
  expect_identical(update(q, parameter = "sigma", remove = 2)$selected, p$selected)
  local_mocked_bindings(prior_fit_candidates = function(...) stop("forced numerical failure"))
  q = p + update_priors("x", belief(mean = 1))
  expect_identical(q$selected, p$selected)
  expect_match(q$candidates$x$message, "mathematical incompatibility was not established", fixed = TRUE)
})

test_that("unsupported brms models and incompatible moment sets fail safely", {
  p = brms_example()
  expect_error(identify_priors(y ~ x + (1 | g), p$model$data, package = "brms"), "Mixed models")
  expect_error(identify_priors(y ~ x, p$model$data, family = gaussian("log"), package = "brms"), "Supported families")
  expect_error(identify_priors(y ~ x + offset(z), p$model$data, package = "brms"), "Offsets")
  q = p + update_priors("sigma", belief(mean = 2), belief(mean = 3))
  expect_identical(q$selected$sigma, p$defaults$sigma)
  expect_identical(q$candidates$sigma$status, "no_supported_fit")
  expect_length(q$beliefs$sigma, 2)
  expect_error(as_brms(q), "allow_fallback")
  q = brms_example(poisson())
  d = q$model$data
  d$y = 2
  expect_s3_class(identify_priors(y ~ x, d, poisson(), package = "brms"), "flexbayes_priors")
})
