# ============================================================
# TEST DATA
# ============================================================

hypothesis_test_draws = function() {
  
  return(cbind(
    advertising = c(0.5, 1.0, 1.5, 2.0),
    price = c(0.8, 0.8, 0.8, 0.8)
  ))
}


# Helper for testing model-dependent functions
evaluate_test_hypothesis = function(expression) {
  
  testthat::local_mocked_bindings(
    bayes_hypothesis_draws = function(object) hypothesis_test_draws(),
    .package = "flexbayes"
  )
  
  return(bayes_evaluate_hypothesis(NULL, expression))
}


# ============================================================
# EXPRESSION VALIDATION
# ============================================================

test_that("valid inequalities are accepted", {
  
  expect_silent(bayes_hypothesis_expression_errors(quote(advertising > 1)))
  expect_silent(bayes_hypothesis_expression_errors(quote(advertising <= 1)))
  expect_silent(bayes_hypothesis_expression_errors(quote(advertising >= price)))
})


test_that("invalid hypotheses are rejected", {
  
  expect_error(bayes_hypothesis_expression_errors(quote(advertising == 1)))
  expect_error(bayes_hypothesis_expression_errors(quote(advertising + 1)))
  expect_error(bayes_hypothesis_expression_errors(quote(advertising)))
  expect_error(bayes_hypothesis_expression_errors(1))
})

test_that("constant coefficient draws are not treated as fixed thresholds", {
  
  draws = hypothesis_test_draws()
  
  result = bayes_hypothesis_plot_parameters(
    draws[, "advertising"],
    draws[, "price"],
    quote(advertising > price)
  )
  
  expect_equal(result$threshold, 0)
  expect_equal(result$plot_draws, draws[, "advertising"] - draws[, "price"])
})

# ============================================================
# EXPRESSION EVALUATION
# ============================================================

test_that("numeric constants are expanded to posterior draws", {
  
  result = bayes_evaluate_expression(1, hypothesis_test_draws())
  
  expect_equal(result, rep(1, 4))
})


test_that("coefficient names return the correct draws", {
  
  result = bayes_evaluate_expression(
    quote(advertising),
    hypothesis_test_draws()
  )
  
  expect_equal(result, c(0.5, 1, 1.5, 2))
})


test_that("binary arithmetic works", {
  
  draws = hypothesis_test_draws()
  
  expect_equal(
    bayes_evaluate_expression(quote(advertising + price), draws),
    draws[, "advertising"] + draws[, "price"]
  )
  
  expect_equal(
    bayes_evaluate_expression(quote(advertising - price), draws),
    draws[, "advertising"] - draws[, "price"]
  )
  
  expect_equal(
    bayes_evaluate_expression(quote(advertising * 2), draws),
    draws[, "advertising"] * 2
  )
  
  expect_equal(
    bayes_evaluate_expression(quote(advertising / 2), draws),
    draws[, "advertising"] / 2
  )
})


test_that("nested arithmetic expressions work", {
  
  draws = hypothesis_test_draws()
  
  result = bayes_evaluate_expression(
    quote((advertising + price) / 2),
    draws
  )
  
  expected = (draws[, "advertising"] + draws[, "price"]) / 2
  
  expect_equal(result, expected)
})


test_that("unary arithmetic works", {
  
  draws = hypothesis_test_draws()
  
  expect_equal(
    bayes_evaluate_expression(quote(-advertising), draws),
    -draws[, "advertising"]
  )
  
  expect_equal(
    bayes_evaluate_expression(quote(+advertising), draws),
    draws[, "advertising"]
  )
})


test_that("invalid expressions produce errors", {
  
  draws = hypothesis_test_draws()
  
  expect_error(
    bayes_evaluate_expression(quote(unknown), draws),
    "Unknown model coefficient"
  )
  
  expect_error(
    bayes_evaluate_expression(quote(log(advertising)), draws),
    "Unsupported operation"
  )
  
  expect_error(
    bayes_evaluate_expression(quote(advertising ^ 2), draws),
    "Unsupported operation"
  )
})


# ============================================================
# PLOTTING PARAMETERS
# ============================================================

test_that("constant thresholds preserve the original scale", {
  
  draws = hypothesis_test_draws()
  left = draws[, "advertising"]
  right = rep(1, nrow(draws))
  
  result = bayes_hypothesis_plot_parameters(
    left, right, quote(advertising > 1)
  )
  
  expect_equal(result$plot_draws, left)
  expect_equal(result$threshold, 1)
  expect_equal(result$x_label, "advertising")
})


test_that("variable thresholds use the difference scale", {
  
  draws = hypothesis_test_draws()
  left = draws[, "advertising"]
  right = draws[, "price"]
  
  result = bayes_hypothesis_plot_parameters(
    left, right, quote(advertising > price)
  )
  
  expect_equal(result$plot_draws, left - right)
  expect_equal(result$threshold, 0)
})


test_that("arithmetic expressions retain meaningful thresholds", {
  
  draws = hypothesis_test_draws()
  
  left = draws[, "advertising"] + draws[, "price"]
  right = rep(2, nrow(draws))
  
  result = bayes_hypothesis_plot_parameters(
    left, right, quote(advertising + price > 2)
  )
  
  expect_equal(result$plot_draws, left)
  expect_equal(result$threshold, 2)
})


# ============================================================
# HYPOTHESIS EVALUATION
# ============================================================

test_that("greater-than hypotheses calculate probabilities", {
  
  result = evaluate_test_hypothesis(quote(advertising > 1))
  
  expect_equal(result$probability, 0.5)
  expect_equal(result$satisfied, c(FALSE, FALSE, TRUE, TRUE))
  expect_s3_class(result, "bayes_hypothesis")
})


test_that("less-than hypotheses calculate probabilities", {
  
  result = evaluate_test_hypothesis(quote(advertising < 1))
  
  expect_equal(result$probability, 0.25)
})


test_that("inclusive inequalities work", {
  
  greater = evaluate_test_hypothesis(quote(advertising >= 1))
  less = evaluate_test_hypothesis(quote(advertising <= 1))
  
  expect_equal(greater$probability, 0.75)
  expect_equal(less$probability, 0.5)
})


test_that("coefficients can be compared", {
  
  result = evaluate_test_hypothesis(quote(advertising > price))
  
  expect_equal(result$probability, 0.75)
})


test_that("arithmetic hypotheses work", {
  
  result = evaluate_test_hypothesis(quote(advertising + price > 2))
  
  expect_equal(result$probability, 0.5)
})


test_that("non-finite values are rejected", {
  
  expect_error(
    evaluate_test_hypothesis(quote(advertising / 0 > 1)),
    "non-finite"
  )
})


# ============================================================
# PLOTTING
# ============================================================

test_that("greater-than hypotheses shade above the threshold", {
  
  result = evaluate_test_hypothesis(quote(advertising > 1))
  data = bayes_hypothesis_plot_data(result)
  
  expect_true(all(data$shaded$x >= 1))
  expect_true(any(data$shaded$x == 1))
})


test_that("less-than hypotheses shade below the threshold", {
  
  result = evaluate_test_hypothesis(quote(advertising < 1))
  data = bayes_hypothesis_plot_data(result)
  
  expect_true(all(data$shaded$x <= 1))
  expect_true(any(data$shaded$x == 1))
})


test_that("plotting returns a ggplot object", {
  
  result = evaluate_test_hypothesis(quote(advertising > 1))
  p = bayes_plot_hypothesis(result)
  
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "advertising")
})


# ============================================================
# PUBLIC FUNCTION
# ============================================================

test_that("hypothesize returns the expected result", {
  
  testthat::local_mocked_bindings(
    bayes_hypothesis_draws = function(object) hypothesis_test_draws(),
    .package = "flexbayes"
  )
  
  result = hypothesize(NULL, advertising > 1)
  
  expect_s3_class(result, "bayes_hypothesis")
  expect_equal(result$probability, 0.5)
})


test_that("hypothesize stores a ggplot when requested", {
  
  testthat::local_mocked_bindings(
    bayes_hypothesis_draws = function(object) hypothesis_test_draws(),
    .package = "flexbayes"
  )
  
  result = hypothesize(NULL, advertising > 1, plot = TRUE)
  
  expect_s3_class(result$plot, "ggplot")
})