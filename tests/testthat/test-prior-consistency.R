test_that("disjoint beliefs are rejected before candidate fitting", {
  p = prior_example()
  p = update(p, parameter = "x", belief(mean = 0),
    belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_identical(p$candidates$x$status, "contradictory")
  expect_equal(p$candidates$x$conflicts, c(2L, 3L))
  expect_null(p$candidates$x$distribution)
  expect_null(p$selected$x)
  expect_length(p$beliefs$x, 3)
  output = paste(capture.output(summary(p, diagnostics = TRUE)), collapse = "\n")
  expect_match(output, "Beliefs 2 and 3 are contradictory", fixed = TRUE)
  expect_match(output, "1.90")
  expect_match(output, "mutually exclusive")
  expect_match(output, "fitting skipped")
  p = update(p, parameter = "x", remove = c(1, 3))
  expect_identical(p$candidates$x$status, "selected")
})

test_that("CDF order and equivalent statements are checked", {
  check = prior_check_probability_constraints
  expect_true(check(list(belief(value < 0, probability = .9), belief(value < 1, probability = .2)))$contradictory)
  expect_false(check(list(belief(value < 0, probability = .2), belief(value < 1, probability = .9)))$contradictory)
  expect_false(check(list(belief(value < 0, probability = .2), belief(value < 0, probability = .2)))$contradictory)
  result = check(list(belief(value < 0, probability = .2), belief(value >= 0, probability = .7)))
  expect_true(result$contradictory)
  expect_match(result$message, "equivalent probability statements disagree")
  expect_false(check(list(belief(value <= 0, probability = .2), belief(value > 0, probability = .8)))$contradictory)
})

test_that("strict inequalities leave room for point masses", {
  check = prior_check_probability_constraints
  expect_false(check(list(belief(value < 0, probability = .2), belief(value <= 0, probability = .8)))$contradictory)
  expect_true(check(list(belief(value <= 0, probability = .8), belief(value < 0, probability = .9)))$contradictory)
  expect_true(check(list(belief(value < 0, probability = .7), belief(value >= 0, probability = .7)))$contradictory)
  expect_false(check(list(belief(value <= 0, probability = .7), belief(value >= 0, probability = .7)))$contradictory)
  expect_false(check(list(belief(value < 0, probability = .4), belief(value > 0, probability = .4)))$contradictory)
  expect_true(check(list(belief(value <= 0, probability = .7), belief(value > 1, probability = .7)))$contradictory)
})

test_that("boundary probabilities are consistent but not normal fits", {
  check = prior_check_probability_constraints
  expect_false(check(list(belief(value < 0, probability = 0), belief(value <= 1, probability = 1)))$contradictory)
  expect_true(check(list(belief(value < 0, probability = 1), belief(value <= 1, probability = 0)))$contradictory)
  expect_false(check(list(belief(value < 0, probability = .5), belief(value < 0, probability = .500001)))$contradictory)
  expect_error(belief(value < 0, probability = -.1))
  p = prior_example()
  p = update(p, parameter = "x", belief(value < 0, probability = 0))
  expect_identical(p$candidates$x$status, "no_supported_fit")
  expect_null(p$selected$x)
  expect_match(p$candidates$x$message, "consistent")
})
