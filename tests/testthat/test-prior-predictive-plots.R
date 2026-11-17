renderer_check = function(family = gaussian(), formula = y ~ x + z + g) {
  prior_example()
  d = data.frame(x = seq(-1, 1, length.out = 36), z = rep(seq(-2, 2, length.out = 6), 6),
    g = factor(rep(c("a", "b", "c"), 12)), y = rep(c(0, 1, 0), 12))
  if (family$family == "gaussian") d$y = 10 + d$x + d$z + rep(c(-1, 1), 18)
  if (family$family == "poisson") d$y = rep(0:3, 9)
  p = identify_priors(formula, d, family)
  for (name in setdiff(names(p$selected), "sigma")) p = update(p, parameter = name, belief(mean = 0), belief(sd = .1))
  set.seed(882)
  check_prior(p, 40)
}

require_draw_renderer = function() {
  skip_if_not("prepare_flexplot" %in% getNamespaceExports("flexplot"), "New public flexplot renderer required")
}

test_that("outcome plots pool stored values with IDs and deterministic bounded display", {
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    before = x
    set.seed(5)
    rng = .Random.seed
    a = plot(x)
    expect_equal(a$data, plot(x, type = "outcome", response = "simulated")$data)
    expect_identical(.Random.seed, rng)
    expect_equal(a$data$value, as.vector(x$simulated))
    expect_equal(a$data$value, x$simulated[cbind(a$data$.draw, a$data$.observation)])
    b = plot(x, response = "expected", max_display = 50)
    expect_equal(nrow(b$data), 50L)
    expect_equal(b$data$value, x$expected[cbind(b$data$.draw, b$data$.observation)])
    expect_identical(x, before)
    expect_s3_class(ggplot2::ggplotGrob(a + ggplot2::theme_bw()), "gtable")
    expect_s3_class(ggplot2::ggplotGrob(b), "gtable")
    expect_s3_class(a$layers[[1]]$stat, if (family$family == "gaussian") "StatBin" else "StatCount")
  }
  x = renderer_check()
  x$simulated = matrix(seq_len(2000000), 1000)
  a = plot(x)
  expect_equal(nrow(a$data), 100000L)
  expect_equal(length(x$simulated), 2000000)
})

test_that("public renderer handles numeric, categorical and conditioned prior draws", {
  require_draw_renderer()
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    for (formula in list(NULL, y ~ x, y ~ g, y ~ x + g, y ~ x | z, y ~ x | g, y ~ g + x | z)) {
      before = x
      p = plot(x, type = "model", formula = formula, draws = 30)
      expect_s3_class(p + ggplot2::theme_bw(), "ggplot")
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
      expect_identical(x, before)
      expect_false(any(vapply(p$layers, function(layer) inherits(layer$stat, "StatSmooth"), logical(1))))
    }
    p = plot(x, type = "model", formula = y ~ x | z, response = "simulated", draws = 4)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    stored = p$data
    expect_equal(stored$y, as.vector(x$simulated[unique(stored$.prior_draw), , drop = FALSE]))
    prep = prior_model_preparation(x, y ~ x | z, 3, NULL, NULL)
    linear = x$coefficient_draws[1:3, , drop = FALSE] %*% t(prior_model_grid_matrix(x, prep$grid))
    expected = switch(family$family, gaussian = linear, binomial = plogis(linear), poisson = exp(linear))
    expect_equal(prior_model_expected(x, prep$grid, 1:3), expected)
  }
})

test_that("new grids preserve contrasts, transformations and coherent coefficient vectors", {
  require_draw_renderer()
  for (formula in list(y ~ x * g + z, y ~ log(x + 2) + g, y ~ poly(x, 2) + g)) {
    x = renderer_check(formula = formula)
    prep = prior_model_preparation(x, y ~ x | g, 3, NULL, NULL)
    ids = c(1, 5, 20)
    values = prior_model_expected(x, prep$grid, ids)
    design = prior_model_grid_matrix(x, prep$grid)
    expect_equal(values, x$coefficient_draws[ids, , drop = FALSE] %*% t(design))
    p = plot(x, type = "model", formula = y ~ x | g, draws = 3)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
  x = renderer_check(formula = y ~ x + g)
  x$metadata$contrasts = list(g = contr.sum(3))
  # A check's design is immutable: incompatible contrast columns must fail.
  expect_error(prior_model_grid_matrix(x, data.frame(x = 0, g = factor("a", levels = c("a", "b", "c")))), "parameterization")
})

test_that("model curves use stored draws, not observed outcomes or new RNG draws", {
  require_draw_renderer()
  x = renderer_check()
  set.seed(52)
  rng = .Random.seed
  p = plot(x, type = "model", formula = y ~ x, draws = 30)
  expect_identical(.Random.seed, rng)
  # Raw observations are hidden; flexplot may preserve their axis extent.
  expect_equal(nrow(p$layers[[1]]$data), 0L)
  altered = x
  altered$metadata$prediction_data$y = 1e50
  q = plot(altered, type = "model", formula = y ~ x, draws = 30)
  lines = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1)))
  expect_equal(ggplot2::ggplot_build(p)$data[[lines[1]]]$y,
    ggplot2::ggplot_build(q)$data[[lines[1]]]$y)
  # Predictor-only newdata is valid: preparation uses stored expected medians.
  priors = prior_example()
  new = check_prior(priors, 10, data.frame(x = c(1, 3, 5), g = c("a", "b", "a")))
  expect_s3_class(ggplot2::ggplotGrob(plot(new, type = "model", formula = y ~ x)), "gtable")
})

test_that("interval aliases and malformed requests remain informative", {
  x = renderer_check()
  expect_equal(plot(x, type = "interval")$data, plot(x, type = "both")$data)
  expect_equal(plot(x, type = "interval", response = "expected")$data, plot(x, type = "expected")$data)
  expect_error(plot(x, max_display = 0), "positive integer")
  expect_error(plot(x, response = "wrong"))
  require_draw_renderer()
  expect_error(plot(x, type = "model", formula = y ~ missing), "Visualization variables")
  expect_error(plot(x, type = "model", formula = y ~ log(x)), "untransformed visualization")
  expect_error(plot(x, type = "model", draws = 0), "positive integer")
  expect_error(plot(x, type = "model", prediction.alpha = 2), "between zero and one")
  expect_error(plot(x, type = "model", response = "simulated", max_display = 2), "Too many prediction cases")
})

test_that("30-draw conditioned smoke plot preserves paths and custom contrasts", {
  require_draw_renderer()
  x = renderer_check()
  p = plot(x, type = "model", formula = y ~ x | z, draws = 30)
  built = ggplot2::ggplot_build(p)
  lines = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1)))
  paths = built$data[[lines[1]]]
  expect_equal(as.integer(tapply(paths$group, paths$PANEL, function(ids) length(unique(ids)))), rep(30L, 3))
  d = data.frame(y = rep(c(0, 1, 0), 12), x = seq(-1, 1, length.out = 36), g = factor(rep(c("a", "b", "c"), 12)))
  contrasts(d$g) = contr.sum(3)
  priors = identify_priors(y ~ x * g, d)
  check = check_prior(priors, 3)
  prep = prior_model_preparation(check, y ~ x | g, 3, NULL, NULL)
  design = model.matrix(delete.response(check$metadata$terms), prep$grid, contrasts.arg = check$metadata$contrasts)
  expect_equal(prior_model_expected(check, prep$grid, 1:3), check$coefficient_draws %*% t(design))
  expect_s3_class(ggplot2::ggplotGrob(plot(check, type = "model", formula = y ~ x | g)), "gtable")
})


test_that("unsupported saved metadata and changed designs fail rather than guess", {
  require_draw_renderer()
  x = renderer_check()
  old = x
  old$metadata$prediction_data = NULL
  expect_error(plot(old, type = "model"), "rerun check_prior")
  x$metadata$prediction_data$x = x$metadata$prediction_data$x + 1
  expect_error(plot(x, type = "model"), "changed since check_prior")
  x = renderer_check(formula = y ~ 1)
  expect_error(plot(x, type = "model"), "predictor")
  expect_s3_class(plot(x), "ggplot")
  x = renderer_check(formula = log(y) ~ x + z)
  expect_error(plot(x, type = "model"), "untransformed response")
})

test_that("outcome labels name the response and count simulations rather than pooled values", {
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    before = x
    x$metadata$formula[[2]] = quote(price)
    p = plot(x, max_display = 50)
    expect_equal(p$labels$title, "Hypothetical distribution of price based on priors")
    expect_equal(p$labels$x, "Simulated price based on priors")
    expect_equal(p$labels$y, if (family$family == "gaussian") "Density" else "Proportion")
    expect_equal(p$labels$caption, "Based on 40 prior draws.")
    expect_null(p$labels$subtitle)
    expect_equal(p$theme, ggplot2::theme_bw())
    q = plot(x, response = "expected", max_display = 20)
    expect_equal(q$labels$title, "Hypothetical distribution of expected price based on priors")
    expect_equal(q$labels$x, "Expected price (based on priors)")
    expect_equal(q$labels$y, "Density")
    expect_equal(q$labels$caption, "Based on 40 prior draws.")
    expect_equal(q$theme, ggplot2::theme_bw())
    expect_s3_class(ggplot2::ggplotGrob(q), "gtable")
    expect_equal(x$expected, before$expected)
    expect_equal(x$simulated, before$simulated)
    expect_equal(x$coefficient_draws, before$coefficient_draws)
  }
  x$metadata$formula[[2]] = as.name("sale price")
  expect_equal(plot(x)$labels$x, "Simulated sale price based on priors")
  x$metadata$formula[[2]] = quote(log(price))
  expect_equal(plot(x)$labels$x, "Simulated log(price) based on priors")
})

test_that("model labels count actually displayed joint draws and apply theme_bw", {
  require_draw_renderer()
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    before = x
    for (count in c(5, 30, 100)) {
      p = plot(x, type = "model", formula = y ~ x | z, draws = count)
      expect_equal(p$labels$title, "Expected relationship based on priors")
      expect_equal(p$labels$y, "y")
      expect_equal(p$labels$caption, paste("Based on", min(count, 40), "prior draws."))
      expect_null(p$labels$subtitle)
      expect_equal(p$theme, ggplot2::theme_bw())
      expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    }
    p = plot(x, type = "model", formula = y ~ x | z, response = "simulated", max_display = 72)
    expect_equal(p$labels$title, "Simulated outcomes based on priors")
    expect_equal(p$labels$caption, "Based on 2 prior draws.")
    expect_equal(p$theme, ggplot2::theme_bw())
    expect_identical(x, before)
  }
})

test_that("interval styling and facet labels preserve original numerical intervals", {
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    before = x
    p = plot(x, type = "interval")
    expect_equal(p$labels$title, "Prior predictive intervals")
    expect_null(p$labels$subtitle)
    expect_equal(p$theme, ggplot2::theme_bw())
    expect_equal(unique(p$data$response), c("Expected outcome", "Predicted outcome"))
    for (kind in c("expected", "simulated")) {
      q = plot(x, type = kind)
      target = t(apply(x[[kind]], 2, quantile, probs = c(.05, .5, .95), names = FALSE))
      expect_equal(as.matrix(q$data[c("lower", "median", "upper")]), target, ignore_attr = TRUE)
    }
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    expect_identical(x, before)
  }
})

test_that("model y axes are independent of observed outcome outliers", {
  require_draw_renderer()
  # Pending a public prediction-only axis-training option in flexplot. Do not
  # remove training layers or replace observed outcomes to make this test pass.
  skip("Blocked: flexplot has no public prediction-only y-axis training control")
  x = renderer_check()
  altered = x
  altered$metadata$prediction_data$y[c(1, 2)] = c(-1e6, 1e6)
  expect_identical(x$coefficient_draws, altered$coefficient_draws)
  expect_identical(x$expected, altered$expected)
  expect_identical(x$metadata$prediction_data[c("x", "z", "g")],
    altered$metadata$prediction_data[c("x", "z", "g")])
  p = ggplot2::ggplot_build(plot(x, type = "model", formula = y ~ x | z))
  q = ggplot2::ggplot_build(plot(altered, type = "model", formula = y ~ x | z))
  expect_equal(lapply(p$layout$panel_params, function(panel) panel$y.range),
    lapply(q$layout$panel_params, function(panel) panel$y.range))
})

test_that("model plots repeat the original observed range across facets and groups", {
  require_draw_renderer()
  for (family in list(gaussian(), binomial(), poisson())) {
    x = renderer_check(family)
    before = x
    for (formula in list(y ~ x | z, y ~ g + x | z, y ~ g | z)) {
      for (response in c("expected", "simulated")) {
        p = plot(x, type = "model", formula = formula, response = response, draws = 4)
        index = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline") && ".observed_value" %in% names(layer$data), logical(1)))
        expect_length(index, 1)
        layer = p$layers[[index]]
        expect_equal(layer$data$.observed_value, range(x$metadata$prediction_data$y))
        expect_equal(layer$aes_params$linetype, "dashed")
        expect_equal(layer$data$.observed_label, c("Observed minimum", "Observed maximum"))
        built = ggplot2::ggplot_build(p)
        references = built$data[[index]]
        expect_equal(length(unique(references$PANEL)), 3L)
        for (panel in unique(references$PANEL)) {
          expect_equal(references$yintercept[references$PANEL == panel], range(x$metadata$observed_response))
        }
        # Removing our reference layers leaves every prediction curve unchanged.
        without = p
        without$layers = without$layers[-c(index, index + 1L)]
        lines = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1)))
        expect_equal(built$data[[lines[1]]]$y, ggplot2::ggplot_build(without)$data[[lines[1]]]$y)
        for (panel in built$layout$panel_params) {
          expect_lte(panel$y.range[1], min(x$metadata$observed_response))
          expect_gte(panel$y.range[2], max(x$metadata$observed_response))
        }
        expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
        expect_equal(p$theme, ggplot2::theme_bw())
        expect_equal(p$labels$caption, "Based on 4 prior draws.")
      }
    }
    expect_identical(x, before)
    for (type in c("outcome", "interval")) {
      p = plot(x, type = type)
      expect_false(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
    }
  }
})

test_that("newdata response values never replace original observed range references", {
  require_draw_renderer()
  priors = prior_example()
  original = model.response(priors$model$frame)
  newdata = data.frame(x = c(1, 3, 5), g = c("a", "b", "a"))
  set.seed(12)
  a = check_prior(priors, 10, newdata)
  newdata$y = c(-1e6, NA, 1e6)
  set.seed(12)
  b = check_prior(priors, 10, newdata)
  expect_identical(a$metadata$observed_response, original)
  expect_identical(b$metadata$observed_response, original)
  expect_identical(a$coefficient_draws, b$coefficient_draws)
  expect_identical(a$expected, b$expected)
  p = plot(a, type = "model", formula = y ~ x)
  index = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline") && ".observed_value" %in% names(layer$data), logical(1)))
  expect_equal(p$layers[[index]]$data$.observed_value, range(original))
  # Finite but irrelevant newdata outcomes can affect flexplot's preparation,
  # but they must never be used as the original observed reference range.
  b$metadata$prediction_data$y = c(-1e6, 0, 1e6)
  q = plot(b, type = "model", formula = y ~ x)
  expect_equal(q$layers[[index]]$data$.observed_value, range(original))
})

test_that("observed references safely handle invalid, constant and legacy observations", {
  require_draw_renderer()
  x = renderer_check()
  x$metadata$observed_response = c(NA, NaN, -Inf, Inf, -3, 7)
  p = plot(x, type = "model", formula = y ~ x)
  index = which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline") && ".observed_value" %in% names(layer$data), logical(1)))
  expect_equal(p$layers[[index]]$data$.observed_value, c(-3, 7))
  x$metadata$observed_response = rep(2, 3)
  p = plot(x, type = "model", formula = y ~ x)
  expect_equal(p$layers[[index]]$data$.observed_value, c(2, 2))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  x$metadata$observed_response = c(NA, Inf)
  expect_warning({ p = plot(x, type = "model", formula = y ~ x) }, "No finite original")
  expect_false(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
  x$metadata$observed_response = NULL
  p = plot(x, type = "model", formula = y ~ x)
  expect_equal(p$layers[[index]]$data$.observed_value, range(x$metadata$prediction_data$y))
  x$metadata$newdata = TRUE
  expect_warning(plot(x, type = "model", formula = y ~ x), "rerun check_prior")
})
