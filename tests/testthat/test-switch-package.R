switch_example = function(package = "rstanarm") {
  skip_if_not_installed("brms")
  p = joint_example()
  if (package == "brms") return(identify_priors(p$model$formula, p$model$data, package = package))
  p
}

test_that("backend switching preserves the model and resolved belief metadata", {
  lower = 2
  p = prior_example() + update_priors("x", belief(mode = 1), belief(sd = 2),
    belief(value < lower, probability = pnorm(.5))) + update_priors("gb", belief(mean = -1))
  skip_if_not_installed("brms")
  original = p
  lower = 100
  set.seed(33)
  seed = .Random.seed
  q = switch_package(p, "brms")
  expect_identical(.Random.seed, seed)
  expect_identical(p, original)
  expect_identical(q$model, p$model)
  expect_identical(q$parameters, p$parameters)
  expect_identical(q$defaults, p$defaults)
  expect_identical(q$default_origin, p$default_origin)
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(q$outcomes, p$outcomes)
  expect_identical(q$next_id, p$next_id)
  expect_identical(q$selected, p$selected)
  expect_equal(q$beliefs$x[[3]]$threshold, 2)
  expect_identical(q$backend$package, "brms")
  expect_identical(q$backend$version, as.character(packageVersion("brms")))
  expect_identical(q$retention_model, prior_retention_model(q))
  expect_type(as_brms(q), "list")
  set.seed(19)
  a = check_prior(p, 50)
  set.seed(19)
  b = check_prior(q, 50)
  expect_identical(a$parameter_draws, b$parameter_draws)
  expect_identical(a$expected, b$expected)
  expect_identical(a$simulated, b$simulated)
  r = switch_package(q, "rstanarm")
  expect_identical(r$selected, p$selected)
  expect_identical(r$beliefs, p$beliefs)
  expect_null(r$backend$coefficients)
  expect_identical(as_rstanarm(r), as_rstanarm(p))
  expect_identical(switch_package(r, "rstanarm"), r)
})

test_that("positive sigma-mode fallbacks recover and incompatible priors cannot leak", {
  # Test package-specific fallback selection independently of automatic migration.
  local_mocked_bindings(prior_maybe_switch_package = identity)
  p = switch_example() + update_priors("sigma", belief(mode = 3), belief(value > 6, probability = .05))
  original = p
  expect_identical(p$candidates$sigma$status, "backend_incompatible")
  q = switch_package(p, "brms")
  expect_identical(q$selected$sigma$distribution, "gamma")
  expect_identical(q$candidates$sigma$status, "selected")
  expect_length(q$fallbacks, 0)
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(q$next_id, p$next_id)
  expect_equal(prior_distribution_operation(q$selected$sigma, "mode"), 3, tolerance = 1e-6)
  expect_s3_class(plot(q, parameter = "sigma"), "ggplot")
  expect_type(as_brms(q), "list")
  expect_s3_class(check_prior(q, 20), "flexbayes_prior_predictive")
  expect_warning({r = switch_package(q, "rstanarm")}, "Switched to rstanarm.*sigma")
  expect_identical(r$selected$sigma, r$defaults$sigma)
  expect_identical(r$fallbacks$sigma$source, "backend_default")
  expect_identical(r$beliefs, q$beliefs)
  expect_error(as_rstanarm(r), "sigma.*allow_fallback")
  expect_error(fit(r), "sigma.*allow_fallback")
  expect_error(check_prior(r), "sigma.*allow_fallback")
  expect_warning({args = as_rstanarm(r, allow_fallback = TRUE)}, "Using retained fallback")
  expect_identical(args$prior_aux, as_rstanarm(switch_example())$prior_aux)
  expect_warning(check_prior(r, 5, allow_fallback = TRUE), "Using retained fallback")
  restored = switch_package(r, "brms")
  expect_identical(restored$selected, q$selected)
  expect_identical(restored$beliefs, q$beliefs)
  expect_identical(p, original)
  removed = restored + update_priors("sigma", remove = 1:2)
  expect_identical(removed$selected$sigma, p$defaults$sigma)
})

test_that("switching enforces whole-model slope families without losing rankings", {
  p = switch_example("brms") + update_priors("x", belief(mean = 0), belief(sd = sqrt(7 / 5)),
    belief(value < 1, probability = pt(1, 7)))
  original = p
  expect_identical(p$selected$x$distribution, "student_t")
  expect_warning({q = switch_package(p, "rstanarm")}, "Switched to rstanarm.*x")
  expect_identical(q$selected$x, q$defaults$x)
  expect_identical(q$selected$z, p$selected$z)
  expect_identical(q$candidates$x$mathematical, p$candidates$x$mathematical)
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(q$candidates$x$status, "backend_incompatible")
  expect_true(q$backend$compatible)
  expect_error(as_rstanarm(q), "allow_fallback")
  expect_identical(switch_package(q, "brms")$selected, p$selected)
  expect_identical(p, original)
  p = p + update_priors("z", belief(mean = 2), belief(sd = .5))
  q = switch_package(p, "rstanarm")
  expect_true(all(vapply(q$selected[c("x", "z")], function(p) p$distribution == "student_t", logical(1))))
  expect_identical(q$beliefs, p$beliefs)
  expect_length(q$fallbacks, 0)
  expect_type(as_rstanarm(q), "list")
})

test_that("outcome reconciliation is recomputed from the same coefficient references", {
  p = switch_example() + update_priors("x", belief(mean = 2), belief(sd = 1)) +
    update_priors(outcome = list(x = 6, z = 1), belief(mean = 12))
  original = p
  q = switch_package(p, "brms")
  expect_identical(q$selected, p$selected)
  expect_identical(q$joint$status, p$joint$status)
  expect_identical(q$joint$outcome, p$joint$outcome)
  expect_identical(q$joint$coefficient, p$joint$coefficient)
  expect_identical(q$outcomes, p$outcomes)
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(q$next_id, p$next_id)
  expect_identical(q$selected$sigma, p$selected$sigma)
  expect_type(as_brms(q), "list")
  expect_identical(switch_package(q, "rstanarm")$selected, p$selected)
  expect_identical(p, original)
  gamma = q + update_priors("sigma", belief(mode = 3))
  expect_warning({back = switch_package(gamma, "rstanarm")}, "sigma.*backend_incompatible")
  expect_identical(back$outcomes, gamma$outcomes)
  expect_identical(back$beliefs, gamma$beliefs)
  expect_identical(back$selected[colnames(back$model$X)], p$selected[colnames(p$model$X)])
  expect_identical(back$selected$sigma, p$defaults$sigma)
  expect_identical(back$joint$status, p$joint$status)
  expect_warning(as_rstanarm(back, allow_fallback = TRUE), "Using retained fallback")
  expect_identical(switch_package(back, "brms")$selected, gamma$selected)
  bad = suppressWarnings(p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 100)))
  q = suppressWarnings(switch_package(bad, "brms"))
  expect_identical(q$joint$status, "incompatible")
  expect_true(q$joint$retained)
  expect_identical(q$selected, bad$selected)
  expect_identical(q$outcomes, bad$outcomes)
  expect_error(as_brms(q), "outcome beliefs.*allow_fallback")
  expect_warning(as_brms(q, allow_fallback = TRUE), "Using retained fallback")
})

test_that("mathematical contradictions and numerical failures remain distinct", {
  p = switch_example() + update_priors("x", belief(mean = 2))
  bad = p + update_priors("x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_warning({q = switch_package(bad, "brms")}, "contradictory")
  expect_identical(q$beliefs, bad$beliefs)
  expect_identical(q$candidates$x$status, "contradictory")
  expect_identical(q$selected$x, p$selected$x)
  expect_identical(q$fallbacks$x$source, "elicited")
  expect_s3_class(plot(q, parameter = "x"), "ggplot")
  expect_error(as_brms(q), "allow_fallback")
  recovered = q + update_priors("x", remove = 2:3)
  expect_length(recovered$fallbacks, 0)
  fit_candidates = prior_fit_candidates
  local_mocked_bindings(prior_fit_candidates = function(...) stop("forced numerical failure"))
  expect_warning({q = switch_package(p, "brms")}, "no_supported_fit")
  expect_match(q$candidates$x$message, "mathematical incompatibility was not established", fixed = TRUE)
  expect_identical(q$selected$x, p$selected$x)
  local_mocked_bindings(prior_fit_candidates = fit_candidates)
  q = switch_package(q, "brms")
  expect_length(q$fallbacks, 0)
  expect_identical(q$selected, p$selected)
})

test_that("GLMs retain link-scale, ratio-scale and standardized constraints", {
  for (family in list(binomial(), poisson())) {
    p = switch_example()
    d = p$model$data
    if (family$family == "binomial") d$y = c(0, 0, 1, 0, 1, 1)
    scale = if (family$family == "binomial") "odds_ratio" else "rate_ratio"
    p = identify_priors(y ~ x + z, d, family) + update_priors("x",
      belief(value < 2, probability = .5, scale = scale), belief(value < 4, probability = .95, scale = scale))
    q = switch_package(p, "brms")
    expect_identical(q$beliefs, p$beliefs)
    expect_identical(q$selected, p$selected)
    expect_type(as_brms(q), "list")
    expect_identical(switch_package(q, "rstanarm")$selected, p$selected)
  }
  p = switch_example() + update_priors("x", belief(mean = 1, scale = "standardized"), belief(sd = .5, scale = "standardized"))
  q = switch_package(p, "brms")
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(q$model$standardization, p$model$standardization)
  expect_identical(q$selected, p$selected)
})

test_that("invalid backend switches reject stale objects atomically", {
  p = switch_example()
  original = p
  for (package in list(NULL, NA_character_, "", "other", 1, c("brms", "rstanarm"))) {
    expect_error(switch_package(p, package), "Choose package")
  }
  expect_error(switch_package(list(), "brms"), "flexbayes_priors")
  altered = p
  altered$model$data$y[1] = 20
  expect_error(switch_package(altered, "brms"), "data changed")
  altered = p
  altered$model$reference[1] = 20
  expect_error(switch_package(altered, "brms"), "metadata changed")
  altered = p
  altered$backend$package = "brms"
  expect_error(switch_package(altered, "rstanarm"), "metadata changed")
  q = p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 12))
  q$outcomes[[1]]$beliefs[[1]]$target = 20
  expect_error(switch_package(q, "brms"), "metadata.*changed")
  q = identify_priors(y ~ x + z, transform(p$model$data, y = 2), poisson(), package = "brms")
  expect_error(switch_package(q, "rstanarm"), "constant responses")
  expect_identical(p, original)
})
