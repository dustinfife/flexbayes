difference_example = function(..., from = list(x = 1, z = 0), to = list(x = 4, z = 0)) {
  update_priors(difference = list(from = from, to = to), ...)
}

difference_moments = function(p, from, to) {
  X = model.matrix(delete.response(terms(p$model$formula)), rbind(as.data.frame(from), as.data.frame(to)))
  contrast = X[2, ] - X[1, ]
  mu = vapply(p$selected[colnames(X)], `[[`, numeric(1), "mean")
  sd = vapply(p$selected[colnames(X)], `[[`, numeric(1), "sd")
  c(mean = sum(contrast * mu), sd = sqrt(sum((contrast * sd)^2)))
}

test_that("differences validate targets, endpoints and beliefs atomically", {
  b = belief(mean = 3)
  expect_error(update_priors(difference = list(), b), "from and to")
  for (d in list(list(from = list(x = 1)), list(from = list(x = 1), to = NULL),
    list(from = list(x = NA), to = list(x = 1)), list(from = list(x = "a"), to = list(x = 1)))) {
    expect_error(update_priors(difference = d, b))
  }
  expect_error(difference_example(remove = integer(), b), "remove applies only")
  expect_error(difference_example(1), "belief")
  expect_error(difference_example(belief(sd = 2)), "mean, mode and probability")
  expect_error(difference_example(belief(mean = 2, scale = "standardized")), "original-scale")
  expect_error(update_priors(parameter = "x", difference = list(), b), "exactly one")
  expect_error(update_priors(outcome = list(x = 1), difference = list(), b), "exactly one")
  expect_error(update_priors("x", difference = list(), b), "exactly one")
  p = joint_example()
  saved = p
  expect_error(p + difference_example(from = list(x = 1), to = list(x = 2), b), "exactly these predictors")
  expect_error(p + difference_example(from = list(x = 1, z = 0), to = list(x = 2, unknown = 1), b), "exactly these predictors")
  expect_identical(p, saved)
  spec = difference_example(belief(mean = 3))
  saved_spec = spec
  q = p + spec
  expect_identical(spec, saved_spec)
  expect_identical(p, saved)
  expect_length(q$outcomes, 1)
})

test_that("the housing contrast implies the exact Normal slope and exports it", {
  skip_if_not_installed("rstanarm")
  from = list(sqft = 1700, age = 44, bedrooms = 3)
  to = list(sqft = 2000, age = 44, bedrooms = 3)
  for (backend in c("rstanarm", "brms")) {
    skip_if_not_installed(backend)
    p = identify_priors(sale_price ~ age + sqft + bedrooms, homes_2024, package = backend)
    q = p + update_priors(difference = list(from = from, to = to),
      belief(mode = 50000), belief(value > 90000, probability = .10))
    expect_identical(q$joint$status, "compatible")
    expect_true(all(q$joint$outcome$satisfied))
    expect_equal(q$selected$sqft$mean, 50000 / 300, tolerance = 1e-6)
    expect_equal(q$selected$sqft$sd, 40000 / qnorm(.9) / 300, tolerance = 1e-5)
    expect_equal(unname(q$joint$design[1, ]), c(0, 0, 300, 0))
    expect_identical(q$selected[c("(Intercept)", "age", "bedrooms", "sigma")],
      p$selected[c("(Intercept)", "age", "bedrooms", "sigma")])
    expect_true(all(vapply(q$joint$stages, prior_joint_stage_verified, logical(1))))
    m = difference_moments(q, from, to)
    expect_lte(abs(pnorm(90000, m[1], m[2], lower.tail = FALSE) - .10), 1e-6)
    if (backend == "rstanarm") {
      exported = as_rstanarm(q)
      expect_equal(unname(exported$prior$location[2]), q$selected$sqft$mean)
      expect_equal(unname(exported$prior$scale[2]), q$selected$sqft$sd)
    } else expect_s3_class(as_brms(q)$prior, "brmsprior")
    set.seed(212)
    check = check_prior(q, n = 1000, newdata = rbind(as.data.frame(from), as.data.frame(to)))
    delta = check$expected[, 2] - check$expected[, 1]
    expect_equal(unname(delta), unname(300 * check$parameter_draws[, "sqft"]), tolerance = 1e-10)
    expect_equal(check$expected, check$coefficient_draws %*% t(check$X))
    text = paste(capture.output(print(summary(q))), collapse = "\n")
    expect_match(text, "Difference in expected sale_price \\(to - from\\)")
    expect_match(text, "From: sqft = 1,700")
    expect_match(text, "P\\(expected difference > 90,000\\)")
    expect_match(paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n"), "Model-matrix contrast")
  }
})

test_that("reversal and changing several predictors preserve shared uncertainty", {
  p = joint_example()
  from = list(x = 1, z = 0)
  to = list(x = 4, z = 1)
  q = p + difference_example(from = from, to = to, belief(mean = 9), belief(value < 12, probability = .9))
  r = p + difference_example(from = to, to = from, belief(mean = -9), belief(value > -12, probability = .9))
  expect_identical(q$joint$status, "compatible")
  expect_equal(q$selected, r$selected, tolerance = 1e-6)
  expect_equal(q$joint$design[1, ], -r$joint$design[1, ])
  expect_equal(unname(q$joint$design[1, ]), c(0, 3, 1))
  m = difference_moments(q, from, to)
  expect_equal(unname(m[2]^2), 9 * q$selected$x$sd^2 + q$selected$z$sd^2)
  expect_identical(q$selected$`(Intercept)`, p$selected$`(Intercept)`)
  expect_identical(q$selected$sigma, p$selected$sigma)
  # Adding endpoint variances would incorrectly include intercept uncertainty.
  X = model.matrix(delete.response(terms(q$model$formula)), rbind(as.data.frame(from), as.data.frame(to)))
  X[, -1] = sweep(X[, -1, drop = FALSE], 2, q$model$reference[-1], "-")
  sd = vapply(q$selected[colnames(X)], `[[`, numeric(1), "sd")
  independent = sum(sweep(X, 2, sd, "*")^2)
  expect_gt(abs(independent - m[2]^2), 1)
})

test_that("differences and outcomes reconcile jointly regardless of update order", {
  p = joint_example()
  a = difference_example(belief(mean = 6))
  b = difference_example(from = list(x = 1, z = 0), to = list(x = 1, z = 1), belief(mode = 2))
  c = update_priors(outcome = list(x = 4, z = 1), belief(mean = 15))
  q = p + a + b + c
  r = p + c + b + a
  expect_identical(q$joint$status, "compatible")
  expect_true(all(q$joint$outcome$satisfied))
  expect_identical(q$selected, r$selected)
  expect_identical(q$joint$baseline, p$selected)
  expect_length(q$outcomes, 3)
  p = p + update_priors("x", belief(mean = 0), belief(sd = 1))
  expect_warning({ q = p + a }, "relaxed")
  expect_identical(q$joint$status, "reconciled")
  expect_identical(q$beliefs, p$beliefs)
  expect_true(all(q$joint$outcome$satisfied))
  expect_false(q$joint$coefficient$satisfied[1])
})

test_that("zero contrasts respect deterministic values and strict comparisons", {
  p = joint_example()
  for (b in list(belief(mean = 0), belief(mode = 0), belief(value < 0, probability = 0),
    belief(value <= 0, probability = 1), belief(value > 0, probability = 0), belief(value >= 0, probability = 1))) {
    q = p + difference_example(from = list(x = 1, z = 0), to = list(x = 1, z = 0), b)
    expect_identical(q$joint$status, "compatible")
    expect_identical(q$selected, p$selected)
    expect_true(all(q$joint$outcome$satisfied))
    expect_silent(prior_joint_validate(q))
  }
  for (b in list(belief(mean = 1), belief(value < 0, probability = 1), belief(value >= 0, probability = 0))) {
    expect_warning({ q = p + difference_example(from = list(x = 1, z = 0), to = list(x = 1, z = 0), b) }, "exactly zero")
    expect_identical(q$joint$status, "incompatible")
    expect_match(q$joint$compatibility$technical$method, "Deterministic")
    expect_identical(q$selected, p$selected)
    expect_error(check_prior(q), "incompatible")
    expect_warning(check_prior(q, n = 2, allow_fallback = TRUE), "fallback")
  }
  q = p + difference_example(from = list(x = 1, z = 0), to = list(x = 1, z = 0), belief(mean = .5, tolerance = 1)) +
    difference_example(belief(mean = 6))
  expect_identical(q$joint$status, "compatible")
  expect_true(all(q$joint$outcome$satisfied))
  expect_warning({ q = p + difference_example(belief(value < 0, probability = .9),
    belief(value < 1, probability = .1)) }, "CDF|probabilit|monotonic")
  expect_identical(q$joint$status, "incompatible")
})

test_that("supported numeric transformations and predictor unit changes use model matrices", {
  p = joint_example(y ~ x * z + I(x^2))
  q = p + difference_example(from = list(x = 1, z = 0), to = list(x = 4, z = 1), belief(mean = 9))
  expect_equal(unname(q$joint$design[1, ]), c(0, 3, 1, 15, 4))
  expect_equal(unname(difference_moments(q, list(x = 1, z = 0), list(x = 4, z = 1))[1]), 9, tolerance = 1e-6)
  p = joint_example()
  q = p + difference_example(belief(mean = 6), belief(value < 9, probability = .9))
  data = p$model$data
  data$x = 1000 * data$x
  r = identify_priors(y ~ x + z, data) + difference_example(from = list(x = 1000, z = 0),
    to = list(x = 4000, z = 0), belief(mean = 6), belief(value < 9, probability = .9))
  expect_equal(q$selected$x$mean, r$selected$x$mean * 1000, tolerance = 1e-6)
  expect_equal(q$selected$x$sd, r$selected$x$sd * 1000, tolerance = 1e-6)
  expect_equal(q$selected$z, r$selected$z, tolerance = 1e-6)
})

test_that("unsupported models fail and numerical failure remains unresolved", {
  p = identify_priors(y ~ x, data.frame(x = 1:6, y = rep(0:1, 3)), family = binomial())
  expect_error(p + difference_example(from = list(x = 1), to = list(x = 2), belief(mean = 1)), "Gaussian identity")
  p = identify_priors(mpg ~ wt + factor(cyl), mtcars)
  expect_error(p + difference_example(from = list(wt = 2, cyl = 4), to = list(wt = 3, cyl = 6), belief(mean = 1)), "numeric predictors")
  p = joint_example(y ~ 0 + x)
  expect_error(p + difference_example(from = list(x = 1), to = list(x = 2), belief(mean = 1)), "intercept")
  p = joint_example()
  local_mocked_bindings(.package = "flexbayes", prior_joint_solve = function(start, ...) list(par = start,
    convergence = 0L, iterations = 1L, violation = 0, stationarity = 1, objective = 0))
  expect_warning({ q = p + difference_example(belief(mean = 6)) }, "could not verify")
  expect_identical(q$joint$status, "unresolved")
  expect_null(q$joint$compatibility)
  expect_identical(q$selected, p$selected)
  expect_error(as_rstanarm(q), "unresolved")
  expect_warning(as_rstanarm(q, allow_fallback = TRUE), "fallback")
})

test_that("multiple probability contrasts and outcomes have order-independent selections", {
  p = joint_example()
  a = difference_example(belief(mode = 6), belief(value < 9, probability = pnorm(1)))
  b = difference_example(from = list(x = 1, z = 0), to = list(x = 1, z = 1),
    belief(mean = 1), belief(value > 2, probability = pnorm(2, lower.tail = FALSE)))
  c = update_priors(outcome = list(x = 3.5, z = .5), belief(mean = 10))
  q = p + a + b + c
  r = p + c + b + a
  expect_identical(q$joint$status, "compatible")
  expect_identical(r$joint$status, "compatible")
  expect_true(all(q$joint$outcome$satisfied))
  expect_identical(q$selected, r$selected)
  expect_identical(q$joint$baseline, p$selected)
  expect_equal(q$selected$x$mean, 2, tolerance = 1e-6)
  expect_equal(q$selected$x$sd, 1, tolerance = 1e-5)
  expect_equal(q$selected$z$mean, 1, tolerance = 1e-6)
  expect_equal(q$selected$z$sd, .5, tolerance = 2e-5)
  set.seed(81)
  saved = .Random.seed
  invisible(p + a)
  expect_identical(.Random.seed, saved)
})

test_that("zero contrasts include transformed coincidences and tolerance boundaries", {
  p = joint_example(y ~ I(x^2))
  q = p + difference_example(from = list(x = -2), to = list(x = 2),
    belief(mean = .5, tolerance = .5), belief(value <= 0, probability = 1))
  expect_identical(q$joint$status, "compatible")
  expect_true(all(q$joint$design == 0))
  expect_identical(q$selected, p$selected)
  expect_true(all(q$joint$outcome$satisfied))
  expect_silent(prior_joint_validate(q))
  expect_warning({ bad = p + difference_example(from = list(x = -2), to = list(x = 2),
    belief(value > 0, probability = .1)) }, "exactly zero")
  expect_identical(bad$joint$status, "incompatible")
  text = paste(capture.output(print(summary(bad, diagnostics = TRUE))), collapse = "\n")
  expect_match(text, "zero contrast")
  expect_match(text, "From: x = -2")
  expect_match(text, "Not verified")
})

test_that("invalid contrast values and unsupported Poisson models fail clearly", {
  for (value in list(Inf, NaN, TRUE, 1 + 1i, c(1, 2), matrix(1), factor("a"))) {
    expect_error(difference_example(from = list(x = value, z = 0), belief(mean = 1)), "finite number")
  }
  for (d in list(list(from = list(x = 1), from = list(x = 2)),
    list(from = list(x = 1), to = list(x = 2), extra = list(x = 3)),
    list(list(x = 1), list(x = 2)))) {
    expect_error(update_priors(difference = d, belief(mean = 1)), "exactly from and to")
  }
  p = identify_priors(y ~ x, data.frame(x = 1:6, y = 1:6), family = poisson())
  expect_error(p + difference_example(from = list(x = 1), to = list(x = 2), belief(mean = 1)), "Gaussian identity")
})

test_that("underflow on a nonzero contrast is not reinterpreted as a zero contrast", {
  q = joint_example() + difference_example(belief(mean = 6))
  problem = q$joint$problem
  n = length(problem$mean)
  theta = c(rep(0, n), rep(-1000, n))
  theta[match("x", names(problem$mean))] = (2 - problem$mean["x"]) / problem$sd["x"]
  expect_equal(unname(prior_joint_values(theta, problem)$value), 6)
  expect_true(any(problem$X != 0))
})
