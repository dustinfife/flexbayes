test_that("brms term extraction unwraps its formula before identifying predictors", {
  skip_if_not_installed("brms")
  model = get_test_model()
  expect_identical(flexplot::get_terms(model), list(predictors = c("hp", "wt"), response = "mpg"))
  original = model
  for (formula in list(mpg ~ hp, mpg ~ hp | wt, mpg ~ wt)) {
    data = bayes_visualize_data(model)
    grid = flexplot:::generate_predictors(data, formula, model)
    expect_setequal(names(grid), c("hp", "wt"))
    omitted = setdiff(c("hp", "wt"), all.vars(formula))
    for (name in omitted) expect_equal(grid[[name]], rep(mean(data[[name]]), nrow(grid)))
    expected = brms::posterior_epred(model, newdata = grid, re_formula = NA)
    expect_equal(bayes_visualize_epred(model, grid), expected)
    predictions = bayes_visualize_predictions(model, data, formula)
    expect_equal(predictions$prediction, colMeans(expected))
    expect_s3_class(flexplot::visualize(model, formula = formula), "ggplot")
  }
  expect_identical(model, original)
})

test_that("omitted factors retain flexplot's representative level and contrasts", {
  skip_if_not_installed("brms")
  model = readRDS(test_path("fixtures", "brms_mtcars_categorical.rds"))
  data = bayes_visualize_data(model)
  grid = flexplot:::generate_predictors(data, mpg ~ hp, model)
  expect_setequal(names(grid), c("hp", "cyl"))
  expected = flexplot:::return_constant_for_predicted_data("cyl", data, model)
  expect_identical(grid$cyl, expected)
  expect_identical(levels(grid$cyl), levels(data$cyl))
  expect_true(all(is.finite(bayes_visualize_epred(model, grid))))
  expect_s3_class(flexplot::visualize(model, formula = mpg ~ hp), "ggplot")
})

test_that("ordinary flexplot term identification is unchanged", {
  model = lm(mpg ~ hp + wt, mtcars)
  expect_identical(flexplot::get_terms(model), list(predictors = c("hp", "wt"), response = "mpg"))
  grid = flexplot:::generate_predictors(mtcars, mpg ~ hp, model)
  expect_setequal(names(grid), c("hp", "wt"))
  expect_equal(grid$wt, rep(mean(mtcars$wt), nrow(grid)))
})

test_that("both Bayesian backends retain all predictors in three-predictor models", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstanarm")
  for (model in list(small_brms, small_stanreg)) {
    original = model
    data = bayes_visualize_data(model)
    expect_identical(flexplot::get_terms(model)$predictors, c("x", "z", "a"))
    for (formula in list(y ~ x, y ~ x + a, y ~ x | z, y ~ x | a)) {
      grid = flexplot:::generate_predictors(data, formula, model)
      expect_setequal(names(grid), c("x", "z", "a"))
      omitted = setdiff(c("x", "z", "a"), all.vars(formula))
      for (name in omitted) expect_identical(grid[[name]],
        flexplot:::return_constant_for_predicted_data(name, data, model))
      predictions = bayes_visualize_predictions(model, data, formula)
      expect_equal(predictions$prediction, unname(colMeans(bayes_visualize_epred(model, grid))))
      plot = flexplot::visualize(model, formula = formula)
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
      plot = flexplot::visualize(model, formula = formula, posterior_lines = 3, seed = 21)
      lines = tail(plot$layers, 1)[[1]]$data
      expect_length(unique(lines$.posterior_draw), 3)
      expect_true(all(is.finite(lines$prediction)))
      rows = match(lines$x, grid$x)
      expect_false(anyNA(rows))
      draws = bayes_visualize_epred(model, grid)
      expect_equal(lines$prediction, unname(draws[cbind(lines$.posterior_draw, rows)]))
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    }
    expect_identical(model, original)
  }
})


test_that("omitted discrete numeric predictors keep flexplot's numeric-mean convention", {
  model = lm(mpg ~ hp + cyl, mtcars)
  grid = flexplot:::generate_predictors(mtcars, mpg ~ hp, model)
  expect_type(grid$cyl, "double")
  expect_equal(grid$cyl, rep(mean(mtcars$cyl), nrow(grid)))
  expect_true(any(grid$cyl != floor(grid$cyl)))
})
