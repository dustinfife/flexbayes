# Tests for Bayesian estimates()
# Uses models saved by data-raw/create-test-models.R


test_that("estimates works with brms and rstanarm", {
  
  data("small_brms", package = "flexbayes")
  data("small_stanreg", package = "flexbayes")
  
  brms_out = flexplot::estimates(small_brms)
  stanreg_out = flexplot::estimates(small_stanreg)
  
  expect_s3_class(brms_out, "bayes_estimates")
  expect_s3_class(stanreg_out, "bayes_estimates")
  
  expect_named(brms_out, c(
    "r.squared", "semi.p", "factor.summary",
    "difference.matrix", "factors", "numbers.summary",
    "numbers", "coefficients", "sigma", "model.comparison"
  ))
  
})


test_that("regression coefficients match posterior draws", {
  
  data("small_brms", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  draws = as.matrix(small_brms)[, "b_x"]
  
  row = out$coefficients[
    out$coefficients$parameter == "x", ,
    drop = FALSE
  ]
  
  expect_equal(row$mean, mean(draws), tolerance = 1e-6)
  expect_equal(row$lower, unname(quantile(draws, 0.025)),
               tolerance = 1e-6)
  
})


test_that("Bayesian R-squared agrees with brms", {
  
  data("small_brms", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  r2 = brms::bayes_R2(small_brms)
  
  expect_equal(
    out$r.squared$mean,
    unname(r2[1, "Estimate"]),
    tolerance = 0.01
  )
  
})


test_that("semi-partial R-squared is reported for each predictor", {
  
  data("small_brms", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  
  expect_setequal(
    out$semi.p$parameter,
    c("x", "z", "a")
  )
  
  expect_true(all(out$semi.p$mean >= 0))
  expect_true(all(out$semi.p$mean <= 1))
  
})


test_that("factor marginal means match posterior predictions", {
  
  data("small_brms", package = "flexbayes")
  data("small_bayes", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  
  newdata = small_bayes
  newdata$a = factor(
    "treatment",
    levels = levels(small_bayes$a)
  )
  
  predictions = brms::posterior_epred(
    small_brms,
    newdata = newdata
  )
  
  expected = mean(rowMeans(predictions))
  
  observed = out$factor.summary$a$mean[
    out$factor.summary$a$parameter == "treatment"
  ]
  
  expect_equal(observed, expected, tolerance = 1e-6)
  
})


test_that("pairwise differences match posterior draws", {
  
  data("small_brms", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  
  # In this additive model, the treatment-control
  # difference equals the treatment coefficient.
  
  draws = as.matrix(small_brms)[, "b_atreatment"]
  
  observed = out$difference.matrix$a$mean[
    out$difference.matrix$a$parameter == "treatment - control"
  ]
  
  expect_equal(observed, mean(draws), tolerance = 1e-6)
  
})


test_that("interaction models are supported", {
  
  data("small_brms_interaction", package = "flexbayes")
  
  out = flexplot::estimates(small_brms_interaction)
  
  expect_true(any(grepl(
    ":",
    out$coefficients$parameter,
    fixed = TRUE
  )))
  
  expect_true("a" %in% names(out$factor.summary))
  
})


test_that("residual standard deviation matches posterior draws", {
  
  data("small_brms", package = "flexbayes")
  
  out = flexplot::estimates(small_brms)
  draws = as.matrix(small_brms)[, "sigma"]
  
  expect_equal(
    out$sigma$mean,
    mean(draws),
    tolerance = 1e-6
  )
  
})