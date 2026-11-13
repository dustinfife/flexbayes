test_that("thresholds resolve caller variables and retain expressions", {
  lower = -2
  b = belief(value < lower, probability = .01)
  expect_equal(b$threshold, lower)
  expect_identical(b$expression, quote(value < lower))
  expect_match(b$label, "value < lower", fixed = TRUE)
  lower = -9
  expect_equal(b$threshold, -2)
  local = local({
    lower = 4
    belief(value >= lower / 2, probability = .5)
  })
  expect_equal(local$threshold, 2)
  expect_identical(local$operator, ">=")
})

test_that("mathematical expressions evaluate only the right hand side", {
  # A promise that would fail if the comparison's left side were evaluated.
  delayedAssign("value", stop("value must not be evaluated"))
  b = belief(value < qt(.01, df = 7), probability = .01)
  expect_equal(b$threshold, qt(.01, df = 7))
  expect_identical(b$expression, quote(value < qt(.01, df = 7)))
  expect_equal(belief(value <= 1 + 2 * 3, probability = .5)$threshold, 7)
  for (operator in c("<", "<=", ">", ">=")) {
    expression = as.call(list(as.name("belief"),
      as.call(list(as.name(operator), as.name("value"), quote(-2.5))), probability = .5))
    b = eval(expression)
    expect_identical(b$operator, operator)
    expect_equal(b$threshold, -2.5)
  }
})

test_that("invalid threshold results and missing variables fail clearly", {
  for (invalid in list(NA_real_, NaN, Inf, -Inf, numeric(), c(1, 2),
      "1", TRUE, NULL, 1 + 2i, matrix(1))) {
    expect_error(belief(value < invalid, probability = .5),
      "single finite numeric value")
  }
  expect_error(belief(value < missing_threshold_variable, probability = .5),
    "Could not evaluate threshold.*missing_threshold_variable")
  expect_error(belief(value < stop("threshold failure"), probability = .5),
    "Could not evaluate threshold: threshold failure", fixed = TRUE)
})

test_that("malformed comparisons fail before evaluating any threshold", {
  expect_error(belief(other < stop("wrong evaluation"), probability = .5), "Use value")
  expect_error(belief(value == stop("wrong evaluation"), probability = .5), "Use value")
  expect_error(belief(value + 1, probability = .5), "Use value")
  expect_error(belief(value, probability = .5), "Use value")
  expect_error(belief((function(a, b) a)(value, 1), probability = .5), "Use value")
  expect_error(belief(value < 1 & value > 0, probability = .5), "Use value")
})
