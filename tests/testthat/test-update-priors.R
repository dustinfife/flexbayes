test_that("specifications are independent of prior objects and validate inputs", {
  spec = update_priors("x", belief(mode = 2), belief(sd = .5))
  expect_s3_class(spec, "flexbayes_prior_update")
  expect_equal(spec$parameter, "x")
  expect_length(spec$beliefs, 2)
  expect_null(spec$beliefs[[1]]$id)
  for (parameter in list(NULL, NA_character_, character(), c("x", "z"), 1, "")) expect_error(update_priors(parameter, belief(mean = 1)), "parameter")
  expect_error(update_priors("x"), "at least one")
  expect_error(update_priors("x", 1), "belief objects")
  expect_error(update_priors("x", belief(sd = 0)), "Invalid belief")
})

test_that("composable addition exactly matches repeated update calls", {
  p = prior_example()
  original = p
  a = update_priors("x", belief(mode = 2), belief(value < 4, probability = .95))
  b = update_priors("gb", belief(mode = -1))
  saved = a
  composed = p + a + b
  repeated = update(p, parameter = "x", belief(mode = 2), belief(value < 4, probability = .95))
  repeated = update(repeated, parameter = "gb", belief(mode = -1))
  expect_identical(composed, repeated)
  expect_identical(p, original)
  expect_identical(a, saved)
  expect_identical(p + a, p + a)
  # Preserve previous beliefs, including their stable IDs.
  appended = composed + update_priors("x", belief(mean = 2))
  expected = update(repeated, parameter = "x", belief(mean = 2))
  expect_identical(appended, expected)
  expect_length(appended$beliefs$x, 3)
  expect_equal(appended$beliefs$x[[1]]$id, composed$beliefs$x[[1]]$id)
})

test_that("existing fitting failures and contradictions retain all constraints", {
  p = prior_example()
  spec = update_priors("x", belief(mode = 150), belief(value < 100, probability = .1), belief(value > 300, probability = .1))
  q = p + spec
  expect_identical(q, update(p, parameter = "x", belief(mode = 150), belief(value < 100, probability = .1), belief(value > 300, probability = .1)))
  expect_length(q$beliefs$x, 3)
  expect_identical(q$candidates$x$status, "no_supported_fit")
  expect_null(q$selected$x)
  spec = update_priors("x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  q = p + spec
  expect_identical(q$candidates$x$status, "contradictory")
  expect_length(q$beliefs$x, 2)
  expect_identical(q, update(p, parameter = "x", belief(value < 0, probability = .95), belief(value > 1, probability = .95)))
})

test_that("addition reuses scale conversion and whole-model compatibility", {
  p = prior_example()
  expect_identical(p + update_priors("x", belief(mode = .5, scale = "standardized")), update(p, parameter = "x", belief(mode = .5, scale = "standardized")))
  q = p + update_priors("x", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  expect_identical(q$candidates$x$status, "backend_incompatible")
  q = q + update_priors("gb", belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
  expect_true(q$backend$compatible)
  expect_equal(q$selected$x$distribution, "student_t")
})

test_that("unsupported operands and unknown parameters fail clearly", {
  p = prior_example()
  spec = update_priors("x", belief(mean = 1))
  expect_error(p + update_priors("missing", belief(mean = 1)), "Unknown parameter")
  expect_error(p + belief(mean = 1), "on the right")
  expect_error(p + 1, "on the right")
  expect_error(1 + spec, "on the left")
  expect_error(spec + p, "on the left")
  expect_error(spec + spec, "on the left")
  expect_error(p + p, "on the right")
  expect_error(+p, "unary")
  expect_error(+spec, "unary")
  # Named beliefs must still be passed as beliefs, not model-control arguments.
  expect_identical(p + update_priors("x", remove = belief(mean = 1)), update(p, parameter = "x", belief(mean = 1)))
})

test_that("thresholds are resolved when specifications are constructed", {
  p = prior_example()
  threshold = 2
  spec = update_priors("x", belief(value < threshold, probability = .95))
  expected = update(p, parameter = "x", belief(value < threshold, probability = .95))
  threshold = 10
  expect_equal(spec$beliefs[[1]]$threshold, 2)
  expect_identical(p + spec, expected)
})
