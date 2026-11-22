grid_example = function(family = gaussian(), package = "rstanarm", n = 25) {
  skip_if_not_installed(package)
  require_draw_renderer()
  d = expand.grid(x = seq(.1, .9, length.out = 13), age = c(10, 20, 30),
    bedrooms = 2:4, g = factor(c("a", "b", "c")))
  d$y = 2 + d$x + sin(seq_len(nrow(d)))
  if (family$family == "binomial") d$y = rep(0:1, length.out = nrow(d))
  if (family$family == "poisson") d$y = rep(0:3, length.out = nrow(d))
  p = identify_priors(y ~ x + age + bedrooms + g, d, family, package)
  for (name in colnames(p$model$X)) p = p + update_priors(name, belief(mean = 0), belief(sd = .03))
  if (family$family == "gaussian") p = p + update_priors("sigma", belief(mean = .2))
  set.seed(938)
  check_prior(p, n)
}


test_that("automatic grid specifications include every predictor without replacing cases", {
  x = grid_example()
  specification = x$metadata$model_grid
  expect_identical(specification$source, "automatic")
  expect_setequal(names(specification$values), c("x", "age", "bedrooms", "g"))
  expect_equal(specification$values$x, seq(.1, .9, length.out = 5))
  expect_equal(specification$values$age, c(10, 20, 30))
  expect_equal(specification$values$bedrooms, 2:4)
  expect_equal(levels(specification$values$g), c("a", "b", "c"))
  prep = prior_model_preparation(x, y ~ x | bedrooms, 3, NULL, NULL)
  expect_equal(length(unique(prep$marginal$data$x)), 51)
  expect_true(all(prep$marginal$data$bedrooms == floor(prep$marginal$data$bedrooms)))
  expect_equal(nrow(prep$marginal$data), 51 * 3 * 3 * 3)
  expect_setequal(prep$marginal$omitted, c("age", "g"))
  expect_equal(ncol(x$expected), nrow(x$metadata$prediction_data))
  expect_equal(plot(x)$data$value, as.vector(x$simulated))
  expect_equal(plot(x, response = "expected")$data$value, as.vector(x$expected))
})

test_that("hypothetical category weights do not follow observed frequencies", {
  require_draw_renderer()
  d = expand.grid(x = seq(.1, .9, length.out = 11), g = factor(c("a", "b", "c")))
  d = rbind(d, d[rep(which(d$g == "a"), 9), ])
  d$y = 1 + d$x + sin(seq_len(nrow(d)))
  p = identify_priors(y ~ x + g, d)
  x = check_prior(p, 3)
  x$coefficient_draws[,] = matrix(rep(c(10, 2, 6, -3), each = 3), 3)
  prep = prior_model_preparation(x, y ~ x, 3, NULL, NULL)
  values = prior_model_expected(x, prep$grid, 1:3)
  expect_equal(as.numeric(values[1, ]), 11 + 2 * prep$grid$x)
  expect_equal(table(attr(prep$grid, "prior_marginal")$cases$g), setNames(rep(51L, 3), c("a", "b", "c")), ignore_attr = TRUE)
  empirical = mean(ifelse(d$g == "b", 6, ifelse(d$g == "c", -3, 0)))
  expect_false(isTRUE(all.equal(as.numeric(values[1, ]), 10 + empirical + 2 * prep$grid$x)))
})

test_that("response-scale marginalization precedes all prior quantiles", {
  for (family in list(gaussian(), binomial(), poisson())) {
    x = grid_example(family, n = 7)
    prep = prior_model_preparation(x, y ~ x | bedrooms, 3, NULL, NULL)
    expected = grid_manual_expected(x, prep)
    actual = prior_model_expected(x, prep$grid, 1:7)
    expect_equal(actual, expected)
    summary = prior_model_summary(x, prep$grid, 20)
    quantiles = t(apply(expected, 2, quantile, c(.05, .5, .95), names = FALSE))
    expect_equal(summary$prediction, unname(quantiles[, 2]))
    expect_equal(summary$.prior_lower, unname(quantiles[, 1]))
    expect_equal(summary$.prior_upper, unname(quantiles[, 3]))
    if (family$family != "gaussian") {
      marginal = attr(prep$grid, "prior_marginal")
      X = prior_model_grid_matrix(x, marginal$cases)
      linear = x$coefficient_draws %*% t(X)
      average = vapply(marginal$groups, function(rows) rowMeans(linear[, rows, drop = FALSE]), numeric(7))
      wrong = if (family$family == "binomial") plogis(average) else exp(average)
      expect_false(isTRUE(all.equal(actual, wrong)))
    }
  }
  d = expand.grid(x = seq(.1, .9, length.out = 11), z = 0:1)
  d$y = 2 + d$x + d$z
  x = check_prior(identify_priors(y ~ x + z, d), 3)
  x$coefficient_draws[,] = rbind(c(0, 0, 3), c(4, 0, -1), c(2, 0, 0))
  prep = prior_model_preparation(x, y ~ x, 3, NULL, NULL)
  summary = prior_model_summary(x, prep$grid, 1000)
  expect_equal(summary$prediction, rep(2, nrow(prep$grid)))
  expect_false(any(summary$prediction == 2.5)) # mean of the two pointwise medians
})

test_that("formula changes and displayed draw counts never resample or alter intervals", {
  x = grid_example(n = 12)
  before = x
  set.seed(21)
  seed = .Random.seed
  for (formula in list(y ~ x, y ~ age | bedrooms, y ~ x | bedrooms, y ~ g + x | age)) {
    a = plot(x, type = "model", formula = formula, display = "both", draws = 2)
    b = plot(x, type = "model", formula = formula, display = "both", draws = 7)
    geom = if (is.numeric(prior_model_preparation(x, formula, 3, NULL, NULL)$plot_data[[all.vars(formula)[2]]])) "GeomRibbon" else "GeomLinerange"
    expect_equal(prior_band_layers(a, geom)[[1]]$data, prior_band_layers(b, geom)[[1]]$data)
    expect_s3_class(ggplot2::ggplotGrob(a), "gtable")
    expect_identical(.Random.seed, seed)
  }
  expect_identical(x, before)
  expect_identical(x$parameter_draws, before$parameter_draws)
})

test_that("explicit newdata is not expanded or replaced by reference profiles", {
  base = grid_example()
  d = base$metadata$prediction_data
  p = identify_priors(y ~ x + age + bedrooms + g, d)
  cases = data.frame(x = c(.2, .2, .8), age = c(10, 30, 20), bedrooms = c(2, 2, 4), g = c("a", "c", "b"))
  x = check_prior(p, 10, cases)
  prep = prior_model_preparation(x, y ~ x | bedrooms, 3, NULL, NULL)
  expect_identical(prep$marginal$source, "explicit")
  expect_equal(nrow(attr(prep$grid, "prior_marginal")$cases), 3)
  expect_equal(nrow(prep$grid), 2)
  expect_equal(prior_model_expected(x, prep$grid, 1:10), grid_manual_expected(x, prep))
  expect_equal(ncol(x$expected), 3)
  expect_equal(plot(x)$data$value, as.vector(x$simulated))
  expect_error(check_prior(p, 10, cases[c("x", "age")]), "missing predictors")
  cases$g[1] = "new"
  expect_error(check_prior(p, 10, cases), "new level")
})

test_that("large grids are bounded, preserve x resolution, and disclose reduction", {
  set.seed(19)
  d = as.data.frame(matrix(runif(50 * 8), 50, 8))
  d$y = rnorm(50)
  p = identify_priors(y ~ ., d)
  x = check_prior(p, 3)
  expect_null(x$metadata$model_grid$data)
  expect_warning({prep = prior_model_preparation(x, y ~ V1, 3, NULL, NULL)}, "reduced numeric resolution")
  expect_lte(nrow(prep$marginal$data), prior_model_grid_policy()$maximum)
  expect_equal(length(unique(prep$marginal$data$V1)), 51)
  expect_true(length(prep$marginal$reduced) > 0)
  expect_equal(as.numeric(prior_model_expected(x, prep$grid, 1)[1, ]), as.numeric(grid_manual_expected(x, prep, 1)))
  d = as.data.frame(matrix(runif(60 * 13), 60, 13))
  d$y = rnorm(60)
  x = check_prior(identify_priors(y ~ ., d), 2)
  expect_error(prior_model_preparation(x, y ~ V1, 3, NULL, NULL), "minimum numeric resolution")
})

test_that("both backends, original factor contrasts and transformations use the grid", {
  for (backend in c("rstanarm", "brms")) {
    x = grid_example(package = backend, n = 5)
    prep = prior_model_preparation(x, y ~ x | bedrooms, 3, NULL, NULL)
    expect_equal(prior_model_expected(x, prep$grid, 1:5), grid_manual_expected(x, prep))
  }
  for (formula in list(y ~ x * g + age, y ~ poly(x, 2) + g + age, y ~ log(x) + g + age)) {
    d = expand.grid(x = seq(.1, .9, length.out = 15), g = factor(c("a", "b", "c")), age = c(10, 20, 30))
    contrasts(d$g) = contr.sum(3)
    d$y = 2 + d$x + sin(seq_len(nrow(d)))
    x = check_prior(identify_priors(formula, d), 5)
    prep = prior_model_preparation(x, y ~ x, 3, NULL, NULL)
    expect_equal(prior_model_expected(x, prep$grid, 1:5), grid_manual_expected(x, prep))
    expect_s3_class(ggplot2::ggplotGrob(plot(x, type = "model", formula = y ~ x)), "gtable")
    expect_equal(prep$marginal$data$g, factor(as.character(prep$marginal$data$g), levels = c("a", "b", "c")))
  }
})

test_that("character levels and invalid joint transformations do not trigger substitution", {
  d = data.frame(x = seq(.1, .9, length.out = 15), g = rep(c("a", "b", "c"), 5))
  d$y = 2 + d$x
  p = identify_priors(y ~ x + g, d)
  x = check_prior(p, 5, newdata = data.frame(x = c(.2, .5), g = c("a", "a")))
  prep = prior_model_preparation(x, y ~ x, 3, NULL, NULL)
  expect_equal(prior_model_expected(x, prep$grid, 1:5), grid_manual_expected(x, prep))
  expect_equal(nrow(prep$marginal$data), 2)
  d = data.frame(x = seq(2, 10, length.out = 20), z = seq(1, 5, length.out = 20))
  d$y = seq_len(20)
  p = identify_priors(y ~ log(x - z), d)
  x = check_prior(p, 5)
  prep = prior_model_preparation(x, y ~ x, 3, NULL, NULL)
  expect_error(suppressWarnings(prior_model_expected(x, prep$grid, 1:5)), "valid configurations")
  explicit = check_prior(p, 5, newdata = d[c(1, 20), c("x", "z")])
  prep = prior_model_preparation(explicit, y ~ x, 3, NULL, NULL)
  expect_equal(prior_model_expected(explicit, prep$grid, 1:5), grid_manual_expected(explicit, prep))
})
