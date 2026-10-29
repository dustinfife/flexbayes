# Bayesian visualization tests

 test_that("brms formulas are handled correctly", {
  
  f = brms::bf(mpg ~ hp + wt)
  
  result = bayes_model_formula(f)
  
  expect_s3_class(result, "formula")
  expect_equal(all.vars(result), c("mpg", "hp", "wt"))
})

test_that("posterior lines respect conditioned panels", {
  
  model = get_test_model()
  
  p = suppressMessages(
    flexplot::visualize(
      model,
      formula = mpg ~ hp | wt,
      posterior_lines = 20,
      seed = 123
    )
  )
  
  expect_length(p$layers, 3)
  
  lines = p$layers[[3]]$data
  
  groups_per_panel = aggregate(
    .posterior_group ~ wt_binned,
    data = lines,
    FUN = function(x) length(unique(x))
  )
  
  expect_equal(nrow(groups_per_panel), 3L)
  expect_true(all(groups_per_panel$.posterior_group == 20))
})

test_that("posterior predictions are numerically correct", {
  
  model = get_test_model()
  data = bayes_visualize_data(model)
  formula = mpg ~ hp
  
  grid = suppressMessages(flexplot:::generate_predictors(
    data = data,
    formula = formula,
    model = model
  ))
  
  expected = colMeans(
    brms::posterior_epred(
      model,
      newdata = grid,
      re_formula = NA
    )
  )
  
  actual = suppressMessages(bayes_visualize_predictions(
    model,
    data,
    formula
  ))
  
  expect_equal(actual$prediction, expected, tolerance = 1e-8)
  expect_equal(actual$model, rep("Posterior mean", nrow(actual)))
})


test_that("visualize produces a valid regression plot", {
  
  model = get_test_model()
  
  p = suppressMessages(flexplot::visualize(
    model,
    plot = "model",
    formula = mpg ~ hp
  ))
  
  expect_s3_class(p, "ggplot")
  
  has_line = any(vapply(
    p$layers,
    function(layer) inherits(layer$geom, "GeomLine"),
    logical(1)
  ))
  
  expect_true(has_line)
})

test_that("default formula generation works", {
  
  model = get_test_model()
  
  p = suppressMessages(flexplot::visualize(
    model,
    plot = "model"
  ))
  
  expect_s3_class(p, "ggplot")
})

test_that("Bayesian plots render correctly across layouts", {
  
  skip_if_not_installed("vdiffr")
  
  model = get_test_model()
  
  plots = list(
    simple = mpg ~ hp,
    multiple = mpg ~ hp + wt,
    conditioned = mpg ~ hp | wt
  )
  
  for (name in names(plots)) {
    
    withr::with_seed(123, {
      
      p = suppressMessages(flexplot::visualize(
        model,
        plot = "model",
        formula = plots[[name]]
      ))
      
      vdiffr::expect_doppelganger(
        paste0("Bayesian ", name),
        p
      )
      
    })
  }
})

test_that("unsupported plot types produce an error", {
  
  model = get_test_model()
  
  expect_error(
    flexplot::visualize(model, plot = "residuals"),
    "only supports plot = 'model'"
  )
})
