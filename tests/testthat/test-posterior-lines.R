test_that("posterior lines preserve draws and conditioned panels", {
  model = get_test_model()
  p = suppressMessages(flexplot::visualize(
    model, formula = mpg ~ hp | wt, posterior_lines = 5, seed = 123
  ))
  expect_s3_class(p, "ggplot")
  lines = p$layers[[length(p$layers)]]$data
  expect_equal(length(unique(lines$.posterior_draw)), 5L)
  expect_true("wt_binned" %in% names(lines))
  expect_gt(length(unique(lines$.posterior_group)), 5L)
})

test_that("posterior line selection is reproducible", {
  model = get_test_model()
  get_lines = function() {
    p = suppressMessages(flexplot::visualize(
      model, formula = mpg ~ hp, posterior_lines = 5, seed = 123
    ))
    p$layers[[length(p$layers)]]$data
  }
  expect_equal(get_lines(), get_lines())
})

test_that("invalid posterior line counts are rejected", {
  model = get_test_model()
  expect_error(flexplot::visualize(model, posterior_lines = -1),
               "nonnegative integer")
  expect_error(flexplot::visualize(model, posterior_lines = 1.5),
               "nonnegative integer")
})

test_that("posterior lines respect grouped layouts", {
  
  model = get_test_model()
  
  p = suppressMessages(
    flexplot::visualize(
      model,
      formula = mpg ~ hp + wt,
      posterior_lines = 20,
      seed = 123
    )
  )
  
  expect_length(p$layers, 3)
  
  lines = p$layers[[3]]$data
  
  expect_true(".posterior_group" %in% names(lines))
  expect_equal(length(unique(lines$.posterior_draw)), 20L)
  
  expect_true(all(is.finite(lines$prediction)))
})

test_that("posterior lines respect categorical groups", {
  
  model = readRDS(
    testthat::test_path(
      "fixtures",
      "brms_mtcars_categorical.rds"
    )
  )
  
  p = suppressMessages(
    flexplot::visualize(
      model,
      formula = mpg ~ hp + cyl,
      posterior_lines = 20,
      seed = 123
    )
  )
  
  expect_length(p$layers, 3)
  
  lines = p$layers[[3]]$data
  
  expect_true("cyl" %in% names(lines))
  expect_equal(length(unique(lines$.posterior_draw)), 20L)
  
  groups_per_draw = aggregate(
    cyl ~ .posterior_draw,
    data = lines,
    FUN = function(x) length(unique(x))
  )
  
  expect_true(all(groups_per_draw$cyl == 3L))
})

test_that("posterior lines render correctly across layouts", {
  
  skip_if_not_installed("vdiffr")
  
  model = get_test_model()
  
  model_categorical = readRDS(
    testthat::test_path(
      "fixtures",
      "brms_mtcars_categorical.rds"
    )
  )
  
  plots = list(
    continuous = list(
      model = model,
      formula = mpg ~ hp
    ),
    conditioned = list(
      model = model,
      formula = mpg ~ hp | wt
    ),
    grouped = list(
      model = model,
      formula = mpg ~ hp + wt
    ),
    categorical = list(
      model = model_categorical,
      formula = mpg ~ hp + cyl
    ),
    categorical_x = list(
      model = model_categorical,
      formula = mpg ~ cyl | hp
    )
  )
  
  for (name in names(plots)) {
    
    withr::with_seed(123, {
      
      p = suppressMessages(
        flexplot::visualize(
          plots[[name]]$model,
          formula = plots[[name]]$formula,
          posterior_lines = 20,
          seed = 123
        )
      )
      
      expect_true(length(p$layers) >= 3L)
      
      posterior_layer = p$layers[[length(p$layers)]]
      
      expect_true(".posterior_group" %in% names(posterior_layer$data))
      expect_equal(
        length(unique(posterior_layer$data$.posterior_draw)),
        20L
      )
      
      vdiffr::expect_doppelganger(
        paste0("Bayesian posterior lines ", name),
        p
      )
    })
  }
})