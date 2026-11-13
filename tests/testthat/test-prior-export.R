test_that("normal export is explicit, ordered, centered and unscaled", {
  p = prior_example()
  p = update(p, parameter = "x", belief(mean = 1.5), belief(sd = .91))
  a = as_rstanarm(p)
  expect_named(a, c("formula", "data", "family", "contrasts", "na.action", "prior", "prior_intercept", "prior_aux", "QR", "sparse"))
  expect_equal(names(a$prior$location), c("x", "gb"))
  expect_equal(unname(a$prior$location), c(1.5, 0))
  expect_equal(unname(a$prior$scale), c(.91, p$selected$gb$sd))
  expect_equal(unname(a$prior_intercept$location), p$selected$`(Intercept)`$mean)
  expect_equal(a$prior_aux$scale, 1 / p$selected$sigma$rate)
  expect_false(a$prior$autoscale)
  expect_false(a$prior_intercept$autoscale)
  expect_false(a$prior_aux$autoscale)
  expect_false(a$QR)
  expect_false(a$sparse)
  expect_equal(a$contrasts, p$model$contrasts)
  expect_equal(a$data, p$model$data)
})

test_that("export matches effective backend preprocessing", {
  p = prior_example()
  captured = new.env()
  testthat::local_mocked_bindings(set_sampling_args = function(..., data) {
    captured$data = data
    stop("captured preprocessing")
  }, .package = "rstanarm")
  expect_error(do.call(rstanarm::stan_glm, as_rstanarm(p)), "captured preprocessing")
  expect_equal(as.numeric(captured$data$prior_scale), c(p$selected$x$sd, p$selected$gb$sd))
  expect_equal(as.numeric(captured$data$xbar), unname(p$model$reference))
  expect_equal(as.numeric(captured$data$prior_scale_for_aux), 1 / p$selected$sigma$rate)
  expect_equal(as.numeric(captured$data$prior_mean_for_intercept), p$selected$`(Intercept)`$mean)
})

test_that("t export uses scale rather than sd and retains df", {
  p = prior_example()
  p = identify_priors(y ~ x, p$model$data)
  p = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  a = as_rstanarm(p)
  expect_equal(a$prior$dist, "t")
  expect_equal(unname(a$prior$df), 7)
  expect_equal(unname(a$prior$scale), 1)
  expect_false(a$prior$autoscale)
  expect_equal(a$prior_intercept$dist, "normal")
  captured = new.env()
  testthat::local_mocked_bindings(set_sampling_args = function(..., data) {
    captured$data = data
    stop("captured t preprocessing")
  }, .package = "rstanarm")
  expect_error(do.call(rstanarm::stan_glm, a), "captured t preprocessing")
  expect_equal(as.numeric(captured$data$prior_df), 7)
  expect_equal(as.numeric(captured$data$prior_scale), 1)
  q = identify_priors(y ~ 0 + x, p$model$data)
  expect_null(as_rstanarm(q)$prior_intercept)
  q = identify_priors(y ~ 1, p$model$data)
  expect_null(as_rstanarm(q)$prior)
})

test_that("fit protects model arguments and passes ordinary arguments unchanged", {
  p = prior_example()
  captured = new.env()
  native = structure(list(), class = c("stanreg", "glm", "lm"))
  testthat::local_mocked_bindings(prior_rstanarm_fit = function(arguments) {
    captured$args = arguments
    native
  })
  expect_identical(fit(p, chains = 2, iter = 600, cores = 1, control = list(adapt_delta = .95)), native)
  expect_equal(captured$args$chains, 2)
  expect_equal(captured$args$control, list(adapt_delta = .95))
  expect_equal(captured$args$prior, as_rstanarm(p)$prior)
  for (name in c("formula", "data", "family", "prior", "prior_aux", "prior_intercept", "QR", "sparse", "contrasts", "weights", "subset", "offset", "na.action", "pri", "unknown")) {
    expect_error(do.call(fit, c(list(p), stats::setNames(list(TRUE), name))), "Cannot override")
  }
  expect_error(fit(p, 2), "unique, complete names")
  expect_error(fit(p, chains = 1, chains = 2), "unique, complete names")
})

test_that("invalid selections and changed models cannot export", {
  p = prior_example()
  q = update(p, parameter = "x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_error(as_rstanarm(q), "contradictory")
  q = update(p, parameter = "x", belief(value < 0, probability = 0))
  expect_error(as_rstanarm(q), "no_supported_fit")
  q = update(p, parameter = "x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  expect_error(as_rstanarm(q), "Incompatible parameters: x.*one slope prior family")
  q = p
  q$model$data$x = q$model$data$x + 1
  expect_error(as_rstanarm(q), "changed since")
  q = p
  q$selected$x$sd = Inf
  expect_error(as_rstanarm(q), "Cannot sample prior")
  q = update(p, parameter = "x", belief(mean = 2))
  q$selected$x$mean = 3
  expect_error(as_rstanarm(q), "does not satisfy recorded beliefs for x")
  q = p
  q$selected$x = NULL
  expect_error(as_rstanarm(q), "missing")
})

test_that("optional Gaussian Stan integration preserves effective priors", {
  skip_if(Sys.getenv("FLEXBAYES_RUN_STAN") != "true", "Set FLEXBAYES_RUN_STAN=true for optional sampling")
  p = prior_example()
  p = identify_priors(y ~ x, p$model$data)
  m = fit(p, chains = 2, iter = 800, warmup = 400, seed = 910, refresh = 0, prior_PD = TRUE)
  expect_s3_class(m, "stanreg")
  draws = as.matrix(m)
  # Backend returns the original intercept: recover its centered prior scale.
  alpha = draws[, "(Intercept)"] + draws[, "x"] * p$model$reference[["x"]]
  expect_equal(mean(alpha), p$selected$`(Intercept)`$mean, tolerance = p$selected$`(Intercept)`$sd * .2)
  expect_equal(sd(draws[, "x"]), p$selected$x$sd, tolerance = p$selected$x$sd * .2)
  expect_equal(mean(draws[, "sigma"]), 1 / p$selected$sigma$rate, tolerance = .2 / p$selected$sigma$rate)
  expect_equal(unname(m$prior.info$prior$scale), unname(p$selected$x$sd))
  expect_null(m$prior.info$prior$adjusted_scale)
  expect_equal(m$prior.info$prior_aux$rate, p$selected$sigma$rate)
  expect_equal(unname(m$prior.info$prior_intercept$scale), p$selected$`(Intercept)`$sd)
  posterior = fit(p, chains = 2, iter = 2000, warmup = 1000, seed = 911, refresh = 0)
  expect_s3_class(posterior, "stanreg")
  expect_equal(unname(posterior$prior.info$prior$scale), p$selected$x$sd)
})

test_that("formula transformations and recorded contrasts survive export", {
  prior_example()
  d = data.frame(y = c(1, 4, 2, 6, 3, 8), x = 1:6, g = factor(rep(c("a", "b"), 3)))
  contrasts(d$g) = stats::contr.sum(2)
  p = identify_priors(y ~ log(x) * g, d)
  a = as_rstanarm(p)
  frame = stats::model.frame(a$formula, a$data)
  expect_equal(stats::model.matrix(a$formula, frame, contrasts.arg = a$contrasts), p$model$X)
  expect_equal(names(a$prior$scale), setdiff(colnames(p$model$X), "(Intercept)"))
  shift = 1
  p = identify_priors(y ~ I(x + shift), d)
  shift = 2
  expect_error(as_rstanarm(p), "changed since")
})
