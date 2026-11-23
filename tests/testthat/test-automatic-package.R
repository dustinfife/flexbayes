automatic_student = function(p, parameter) {
  p + update_priors(parameter, belief(mean = 0), belief(sd = sqrt(7 / 5)),
    belief(value < 1, probability = pt(1, 7)))
}

automatic_sigma = function(p = joint_example()) {
  skip_if_not_installed("brms")
  p + update_priors("sigma", belief(mode = 3))
}

test_that("positive sigma modes migrate once and retain every original answer", {
  p = joint_example() + update_priors("x", belief(mean = 2, tolerance = .001), belief(sd = 1))
  before = p
  set.seed(29)
  seed = .Random.seed
  expect_message({q = automatic_sigma(p)}, "rstanarm.*sigma.*Switching to brms")
  expect_identical(.Random.seed, seed)
  expect_identical(p, before)
  expect_identical(q$backend$package, "brms")
  expect_identical(q$model, p$model)
  expect_identical(q$defaults, p$defaults)
  expect_identical(q$default_origin, p$default_origin)
  expect_identical(q$beliefs$x, p$beliefs$x)
  expect_identical(q$beliefs$sigma[[1]]$id, p$next_id)
  expect_identical(q$next_id, p$next_id + 1L)
  expect_identical(q$candidates$sigma$status, "selected")
  expect_identical(q$selected$sigma$distribution, "gamma")
  expect_equal(prior_distribution_operation(q$selected$sigma, "mode"), 3, tolerance = 1e-6)
  expect_length(q$fallbacks, 0)
  expect_length(q$package_switches, 1)
  expect_identical(q$package_switches[[1]]$parameters, "sigma")
  expect_message(capture.output(print(summary(q))), NA)
  ordinary = paste(capture.output(print(summary(q))), collapse = "\n")
  expect_match(ordinary, "brms")
  expect_false(grepl("Automatic package switch", ordinary, fixed = TRUE))
  expect_match(paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n"), "Automatic package switch: rstanarm -> brms", fixed = TRUE)
  expect_message({r = q + update_priors("sigma", remove = 1)}, NA)
  expect_identical(r$backend$package, "brms")
  expect_identical(r$selected$sigma, p$defaults$sigma)
  expect_identical(r$package_switches, q$package_switches)
})

test_that("migration recovers mixed coefficient families without changing defaults", {
  skip_if_not_installed("brms")
  p = prior_example() + update_priors("gb", belief(mode = 2), belief(sd = .4),
    belief(value < 2.4, probability = pnorm(1)))
  expect_message({q = automatic_student(p, "x")}, "rstanarm.*x.*Switching to brms")
  expect_identical(q$selected$x$distribution, "student_t")
  expect_identical(q$selected$gb, p$selected$gb)
  expect_identical(q$selected$sigma, p$selected$sigma)
  expect_identical(q$beliefs$gb, p$beliefs$gb)
  expect_identical(q$candidates$x$mathematical, q$selected$x)
  expect_type(as_brms(q), "list")
})

test_that("complete outcome reconciliation is rechecked against the same baseline", {
  p = joint_example() + update_priors("x", belief(mean = 2), belief(sd = 1)) +
    update_priors(outcome = list(x = 6, z = 1), belief(mean = 12))
  expect_message({q = automatic_sigma(p)}, "Switching to brms")
  expect_identical(q$outcomes, p$outcomes)
  expect_identical(q$beliefs$x, p$beliefs$x)
  expect_identical(q$joint$outcome, p$joint$outcome)
  expect_identical(q$joint$coefficient, p$joint$coefficient)
  expect_identical(q$selected[colnames(p$model$X)], p$selected[colnames(p$model$X)])
  expect_true(q$joint$status %in% c("compatible", "reconciled"))
  expect_false(isTRUE(q$joint$retained))
  expect_type(as_brms(q), "list")
})

test_that("ordinary and contradictory updates do not attempt a package switch", {
  p = joint_example()
  local_mocked_bindings(prior_package_available = function(...) stop("Unexpected package trial"))
  expect_message({q = p + update_priors("x", belief(mode = 2))}, NA)
  expect_identical(q$backend$package, "rstanarm")
  q = p + update_priors("sigma", belief(value < 0, probability = .5))
  expect_identical(q$candidates$sigma$status, "contradictory")
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$selected, p$selected)
  expect_error(check_prior(q), "contradictory.*allow_fallback")
  local_mocked_bindings(prior_fit_candidates = function(...) stop("forced numerical failure"))
  q = p + update_priors("x", belief(mean = 1))
  expect_identical(q$candidates$x$status, "no_supported_fit")
  expect_match(q$candidates$x$message, "mathematical incompatibility was not established", fixed = TRUE)
  expect_identical(q$backend$package, "rstanarm")
})

test_that("unavailable brms retains priors and requires explicit fallback approval", {
  p = joint_example()
  local_mocked_bindings(prior_package_available = function(package) FALSE)
  expect_message({q = p + update_priors("sigma", belief(mode = 3))}, "rstanarm.*sigma.*Install brms.*evaluate an alternative")
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$selected, p$selected)
  expect_null(q$package_switches)
  expect_identical(q$fallbacks$sigma$status, "backend_incompatible")
  expect_error(check_prior(q), "allow_fallback")
  expect_error(as_rstanarm(q), "allow_fallback")
  expect_warning(check_prior(q, 3, allow_fallback = TRUE), "Using retained fallback")
  expect_warning(as_rstanarm(q, allow_fallback = TRUE), "Using retained fallback")
})

test_that("unverified destination fits and unexportable priors cannot be accepted", {
  skip_if_not_installed("brms")
  p = joint_example()
  raw = prior_apply_parameter(p, belief(mode = 3), parameter = "sigma")
  local_mocked_bindings(as_brms = function(...) stop("unsupported export"))
  expect_message({q = p + update_priors("sigma", belief(mode = 3))}, NA)
  expect_identical(q, raw)
  expect_identical(p$beliefs$sigma, list())
})

test_that("destination numerical failures preserve the original package fallback", {
  skip_if_not_installed("brms")
  p = joint_example()
  raw = prior_apply_parameter(p, belief(mode = 3), parameter = "sigma")
  original = prior_fit_candidates
  local_mocked_bindings(prior_fit_candidates = function(records, reference, backend) {
    if (backend == "brms") stop("forced destination search failure")
    original(records, reference, backend)
  })
  expect_message({q = p + update_priors("sigma", belief(mode = 3))}, NA)
  expect_identical(q, raw)
  expect_identical(q$candidates$sigma$status, "backend_incompatible")
})

test_that("neither package satisfying all sigma beliefs leaves the prior unchanged", {
  skip_if_not_installed("brms")
  p = joint_example()
  specification = update_priors("sigma", belief(mode = 3), belief(mean = 1), belief(sd = .01))
  raw = do.call(prior_apply_parameter, c(list(object = p, parameter = "sigma"), specification$beliefs))
  expect_message({q = p + specification}, NA)
  expect_identical(q, raw)
  expect_length(q$beliefs$sigma, 3)
  expect_error(check_prior(q), "allow_fallback")
})

test_that("unsuccessful outcomes prohibit migration even with a package restriction", {
  skip_if_not_installed("brms")
  p = joint_example() + update_priors(outcome = list(x = 6, z = 1), belief(mean = 12))
  expect_message({q = automatic_student(p, "x")}, NA)
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$candidates$x$status, "backend_incompatible")
  expect_identical(q$joint$status, "unresolved")
  expect_identical(q$selected, p$selected)
  expect_identical(q$outcomes, p$outcomes)
  bad = suppressWarnings(p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 100)))
  expect_message(expect_warning({q = bad + update_priors("sigma", belief(mode = 3))}, "incompatible"), NA)
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$joint$status, "incompatible")
  expect_identical(q$selected, bad$selected)
  expect_error(as_rstanarm(q), "allow_fallback")
})

test_that("explicit switches honor the requested package without automatic reversal", {
  expect_false("switch_backend" %in% getNamespaceExports("flexbayes"))
  expect_true("switch_package" %in% getNamespaceExports("flexbayes"))
  expect_message({p = automatic_sigma()}, "Switching to brms")
  before = p
  expect_message(expect_warning({q = switch_package(p, "rstanarm")}, "sigma.*backend_incompatible"), NA)
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$beliefs, p$beliefs)
  expect_identical(p, before)
  expect_error(check_prior(q), "allow_fallback")
  expect_message({r = q + update_priors("x", belief(mode = 2))}, "Switching to brms")
  expect_identical(r$backend$package, "brms")
})

test_that("automatic selections agree with simulation, export and fit dispatch", {
  expect_message({p = automatic_sigma()}, "Switching to brms")
  set.seed(931)
  expect_message({x = check_prior(p, 30)}, NA)
  expect_identical(x$metadata$backend, p$backend)
  expect_true(all(x$parameter_draws[, "sigma"] > 0))
  set.seed(931)
  for (name in names(p$selected)) {
    prior = p$selected[[name]]
    values = if (name == "sigma") rgamma(30, shape = prior$shape, rate = prior$rate) else
      rnorm(30, prior$mean, prior$sd)
    expect_equal(x$parameter_draws[, name], values)
  }
  expect_message({arguments = as_brms(p)}, NA)
  expect_match(arguments$prior$prior[arguments$prior$class == "sigma"], "gamma", fixed = TRUE)
  expect_s3_class(plot(p, parameter = "sigma"), "ggplot")
  expect_s3_class(plot(x, type = "model"), "ggplot")
  captured = new.env()
  local_mocked_bindings(prior_brms_fit = function(arguments) { captured$arguments = arguments; "sentinel" })
  expect_message(expect_identical(fit(p, chains = 1, sample_prior = "only"), "sentinel"), NA)
  expect_identical(captured$arguments$prior, arguments$prior)
  expect_identical(captured$arguments$formula, p$model$formula)
})


test_that("an acceptable shared-family candidate avoids unnecessary migration", {
  p = prior_example() + update_priors("gb", belief(mode = 2), belief(sd = .4))
  expect_message({q = automatic_student(p, "x")}, NA)
  expect_identical(q$backend$package, "rstanarm")
  expect_identical(q$selected$x$distribution, "student_t")
  expect_identical(q$selected$gb$distribution, "student_t")
  expect_identical(q$beliefs$gb, p$beliefs$gb)
  expect_length(q$fallbacks, 0)
})

test_that("migration preserves converted scales and resolved threshold expressions", {
  skip_if_not_installed("brms")
  p = joint_example() + update_priors("x", belief(mode = .5, scale = "standardized", tolerance = .001), belief(sd = 1))
  upper = prior_distribution_operation(p$selected$x, "quantile", .95)
  p = p + update_priors("x", belief(value < upper, probability = .95))
  records = p$beliefs
  upper = 1e8
  expect_message({q = automatic_sigma(p)}, "Switching to brms")
  expect_identical(q$beliefs$x, records$x)
  expect_identical(q$model$standardization, p$model$standardization)
  expect_identical(q$beliefs$x[[1]]$scale, "standardized")
  expect_identical(q$beliefs$x[[1]]$tolerance, p$beliefs$x[[1]]$tolerance)
  expect_identical(q$beliefs$x[[3]]$threshold, p$beliefs$x[[3]]$threshold)
  expect_identical(q$beliefs$x[[3]]$expression, p$beliefs$x[[3]]$expression)
})

test_that("outcome additions can resolve a package-specific specification", {
  expect_message({p = automatic_sigma()}, "Switching to brms")
  expect_warning({p = switch_package(p, "rstanarm")}, "backend_incompatible")
  records = p$beliefs
  expect_message({q = p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 12))}, "Switching to brms")
  expect_identical(q$backend$package, "brms")
  expect_identical(q$beliefs, records)
  expect_true(q$joint$status %in% c("compatible", "reconciled"))
  expect_length(q$fallbacks, 0)
  expect_type(as_brms(q), "list")
})
