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
