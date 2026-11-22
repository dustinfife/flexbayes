test_that("exactly one target is required and control names match exactly", {
  b = belief(mode = 2)
  expect_error(update_priors(), "exactly one of parameter or outcome")
  expect_error(update_priors(b), "exactly one of parameter or outcome")
  expect_error(update_priors(parameter = "x", outcome = list(x = 1), b), "exactly one")
  expect_error(update_priors("x", outcome = list(x = 1), b), "exactly one")
  expect_error(update_priors("x", parameter = "x", b), "parameter only once")
  expect_error(update_priors(param = "x", b), "exact argument names")
  expect_error(update_priors(out = list(x = 1), b), "exact argument names")
  expect_error(update_priors(parameter = NULL, b), "parameter must")
  expect_error(update_priors(outcome = NULL, b), "named list")
  expect_error(update_priors("x", b, rem = 1), "belief objects")
  for (remove in list(integer(), NULL, 1, NA)) {
    expect_error(update_priors(outcome = list(x = 1), b, remove = remove), "remove applies only to parameter")
  }
  positional = update_priors("x", b, belief(sd = 1))
  named = update_priors(parameter = "x", b, belief(sd = 1))
  expect_identical(positional, named)
  expected = update_priors(outcome = list(x = 1), b, belief(value < 3, probability = .95))
  expect_identical(expected$target, "outcome")
  expect_identical(expected$at, list(x = 1))
  expect_identical(expected$beliefs[[1]], b)
  expect_length(expected$beliefs, 2)
  expect_identical(class(positional), class(expected))
  expect_identical(update_priors(b, outcome = list(x = 1)), update_priors(outcome = list(x = 1), b))
})

test_that("outcome validation and atomic application preserve the shared restrictions", {
  for (outcome in list(1, list(1), list(x = 1, x = 2), setNames(list(1), ""))) {
    expect_error(update_priors(outcome = outcome, belief(mean = 1)), "named list")
  }
  for (value in list(NA, Inf, NaN, "1", TRUE, 1 + 1i, c(1, 2), matrix(1))) {
    expect_error(update_priors(outcome = list(x = value), belief(mean = 1)), "finite number")
  }
  expect_error(update_priors(outcome = list(x = 1)), "belief")
  expect_error(update_priors(outcome = list(x = 1), 1), "belief")
  expect_error(update_priors(outcome = list(x = 1), belief(sd = 1)), "mean, mode and probability")
  expect_error(update_priors(outcome = list(x = 1), belief(mean = 1, scale = "standardized")), "original-scale")
  p = joint_example()
  original = p
  spec = update_priors(outcome = list(x = 1), belief(mean = 1))
  saved = spec
  expect_error(p + spec, "exactly these predictors")
  expect_identical(p, original)
  expect_identical(spec, saved)
  expect_error(p + update_priors(parameter = "missing", belief(mean = 1)), "Unknown parameter")
  expect_error(p + update_priors(parameter = "x", remove = 1), "Invalid removal")
  expect_identical(p, original)
  expect_error(1 + spec, "on the left")
  expect_error(spec + p, "on the left")
  expect_error(spec + spec, "on the left")
  expect_error(+spec, "unary")
  invalid = spec
  invalid$target = "other"
  expect_error(p + invalid, "Unknown update target")
})

test_that("parameter and outcome targets compose through the original mathematical paths", {
  for (backend in c("rstanarm", "brms")) {
    skip_if_not_installed(backend)
    p = joint_example()
    if (backend == "brms") p = identify_priors(p$model$formula, p$model$data, package = backend)
    original = p
    q = p + update_priors(parameter = "x", belief(mode = 2), belief(sd = 1))
    expect_identical(q, prior_apply_parameter(p, parameter = "x", belief(mode = 2), belief(sd = 1)))
    at = list(x = 3, z = 0)
    mu = q$selected$`(Intercept)`$mean + 2 * (at$x - q$model$reference["x"])
    first = update_priors(outcome = at, belief(mode = unname(mu)))
    second = update_priors(outcome = list(x = 1, z = 0), belief(mean = unname(mu - 4)))
    a = q + first + second
    expected = prior_apply_outcome(prior_apply_outcome(q, first), second)
    expect_identical(a, expected)
    expect_identical(a$selected, q$selected)
    expect_identical(a$selected$sigma, p$selected$sigma)
    expect_length(a$outcomes, 2)
    expect_identical(names(a$outcomes[[1]]), c("at", "beliefs"))
    expect_identical(a$outcomes[[1]]$beliefs[[1]]$id, 3L)
    expect_identical(a$outcomes[[2]]$beliefs[[1]]$id, 4L)
    expect_null(first$beliefs[[1]]$id)
    expect_identical(p, original)
    revised = a + update_priors("x", belief(mean = 2), remove = 1)
    expect_identical(revised, prior_apply_parameter(a, parameter = "x", belief(mean = 2), remove = 1))
    expect_identical(revised$outcomes, a$outcomes)
    expect_equal(vapply(revised$beliefs$x, `[[`, numeric(1), "id"), c(2, 5))
    export = if (backend == "brms") as_brms else as_rstanarm
    expect_type(export(revised), "list")
    set.seed(27)
    check = check_prior(revised, 10)
    expect_equal(dim(check$expected), c(10, nrow(p$model$X)))
  }
  intercept = joint_example(y ~ 1) + update_priors(outcome = list(), belief(mean = 2))
  expect_true(all(intercept$joint$outcome$satisfied))
})

test_that("obsolete public belief-update names and methods are removed", {
  ns = asNamespace("flexbayes")
  expect_false("update_outcome" %in% getNamespaceExports("flexbayes"))
  expect_false(exists("update_outcome", ns, inherits = FALSE))
  expect_false(exists("update.flexbayes_priors", ns, inherits = FALSE))
  expect_null(getS3method("+", "flexbayes_outcome_update", optional = TRUE))
  expect_null(getS3method("update", "flexbayes_priors", optional = TRUE))
  expect_false("prior_apply_parameter" %in% getNamespaceExports("flexbayes"))
  expect_false("prior_outcome_update" %in% getNamespaceExports("flexbayes"))
})
