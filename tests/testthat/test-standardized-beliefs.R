test_that("standardized moments convert to original units and preserve inputs", {
  p = prior_example()
  factor = sd(p$model$data$y) / sd(p$model$data$x)
  for (type in c("mean", "mode", "sd")) {
    b = do.call(belief, c(stats::setNames(list(.5), type), list(scale = "standardized")))
    q = update(p, parameter = "x", b)
    r = q$beliefs$x[[1]]
    expect_equal(r$target, .5 * factor)
    expect_equal(r$entered$target, .5)
    expect_equal(r$tolerance, 1e-6 * factor)
    expect_equal(r$standardization_factor, factor)
    expect_equal(prior_distribution_operation(q$selected$x, type), .5 * factor)
  }
  expect_error(belief(mean = 1, scale = "unknown"), "scale")
})

test_that("mixed beliefs and probabilities fit on original units", {
  p = prior_example()
  factor = p$model$standardization$x$factor
  q = update(p, parameter = "x", belief(mode = .5, scale = "standardized"),
    belief(value < 1, probability = .95, scale = "standardized"))
  expect_equal(q$selected$x$mean, .5 * factor)
  expect_equal(q$selected$x$sd, .5 * factor / qnorm(.95))
  expect_equal(q$beliefs$x[[2]]$entered$threshold, 1)
  expect_equal(q$beliefs$x[[2]]$threshold, factor)
  expect_equal(q$beliefs$x[[2]]$tolerance, 1e-6)
  expect_identical(q$beliefs$x[[2]]$expression, quote(value < 1))
  mixed = update(p, parameter = "x", belief(mean = .5, scale = "standardized"), belief(sd = 2))
  expect_equal(mixed$selected$x$mean, .5 * factor)
  expect_equal(mixed$selected$x$sd, 2)
  conflict = update(p, parameter = "x", belief(value < 1, probability = .95, scale = "standardized"), belief(value > 2 * factor, probability = .95))
  expect_identical(conflict$candidates$x$status, "contradictory")
  q = update(q, parameter = "x", remove = 1)
  expect_equal(q$beliefs$x[[1]]$id, 2L)
  expect_equal(q$beliefs$x[[1]]$entered$threshold, 1)
  expect_equal(q$beliefs$x[[1]]$threshold, factor)
})

test_that("summary, standardized densities, export and fit agree", {
  p = prior_example()
  q = update(p, parameter = "x", belief(mean = .5, scale = "standardized"), belief(sd = .2, scale = "standardized"))
  factor = q$model$standardization$x$factor
  text = paste(capture.output(summary(q, diagnostics = TRUE)), collapse = "\n")
  expect_match(text, "mean = 0.5 [standardized]", fixed = TRUE)
  expect_match(text, "SD(y)/SD(x)", fixed = TRUE)
  expect_match(text, "original-scale equivalent")
  plot = plot(q, parameter = "x", scale = "standardized")
  normal = plot$data[grepl("normal", plot$data$Candidate), ]
  expect_equal(normal$density, dnorm(normal$value, .5, .2))
  expect_s3_class(plot + ggplot2::theme_bw(), "ggplot")
  expect_equal(unname(as_rstanarm(q)$prior$location[1]), .5 * factor)
  expect_equal(unname(as_rstanarm(q)$prior$scale[1]), .2 * factor)
  testthat::local_mocked_bindings(prior_rstanarm_fit = function(arguments) arguments)
  arguments = fit(q, iter = 100)
  expect_equal(arguments$data, p$model$data)
  expect_equal(arguments$formula, p$model$formula)
  expect_equal(arguments$prior$scale, as_rstanarm(q)$prior$scale)
  set.seed(12)
  a = check_prior(q, 3)
  expect_equal(a$metadata$priors$x, q$selected$x)
  expect_equal(a$metadata$reference, p$model$reference)
})

test_that("unsupported terms and changed data are rejected", {
  p = prior_example()
  b = belief(mean = .5, scale = "standardized")
  for (name in c("(Intercept)", "gb", "sigma")) expect_error(update(p, parameter = name, b))
  d = p$model$data
  for (formula in list(y ~ log(x), y ~ poly(x, 2), y ~ x * g, y ~ I(x^2))) {
    q = identify_priors(formula, d)
    name = setdiff(names(q$selected), c("(Intercept)", "sigma"))[1]
    expect_error(update(q, parameter = name, b), "untransformed numeric slope")
  }
  q = update(p, parameter = "x", b)
  q$model$data$x = q$model$data$x + 1
  expect_error(update(q, parameter = "x", belief(sd = 1)), "changed since")
  expect_error(check_prior(q), "changed since")
  expect_error(as_rstanarm(q), "changed since")
  expect_error(plot(q, parameter = "x", scale = "standardized"), "changed since")
  expect_error(plot(p, scale = "standardized"), "Select one")
  d$x = 1
  expect_error(identify_priors(y ~ x, d), "full-rank")
  q = identify_priors(y ~ 0 + x, d)
  expect_error(update(q, parameter = "x", b), "positive finite")
  d$y = 1
  expect_error(identify_priors(y ~ 1, d), "positive sd")
  d$y[1] = NA
  expect_error(identify_priors(y ~ 1, d), "missing")
})

test_that("standardization preserves Student t scale and reusable beliefs", {
  base = prior_example()
  p = identify_priors(y ~ x, base$model$data)
  factor = p$model$standardization$x$factor
  p = update(p, parameter = "x", belief(mean = 0, scale = "standardized"),
    belief(sd = sqrt(7 / 5), scale = "standardized"),
    belief(value <= 1, probability = pt(1, 7), scale = "standardized"))
  expect_equal(p$selected$x$distribution, "student_t")
  expect_equal(p$selected$x$scale, factor)
  expect_equal(unname(as_rstanarm(p)$prior$scale), factor)
  density = plot(p, parameter = "x", scale = "standardized")$data
  expect_equal(density$density, dt(density$value, 7))
  # Reusing a stored belief must not convert it twice.
  q = update(p, parameter = "x", p$beliefs$x[[2]])
  expect_equal(q$beliefs$x[[4]]$target, p$beliefs$x[[2]]$target)
  expect_equal(q$beliefs$x[[4]]$entered$target, sqrt(7 / 5))
})
