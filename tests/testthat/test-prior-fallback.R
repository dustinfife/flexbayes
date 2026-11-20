fallback_student_beliefs = function(object, parameter) {
  update(object, parameter = parameter, belief(mean = 0), belief(sd = sqrt(7 / 5)),
    belief(value < 1, probability = pt(1, 7)))
}

test_that("initial unsuccessful sigma elicitation retains the backend default", {
  p = identify_priors(sale_price ~ sqft + age + bedrooms, homes_2024)
  before = p
  q = p + update_priors("sigma", belief(value < 30000, probability = .5), belief(value < 60000, probability = .95))
  expect_identical(q$selected, p$selected)
  expect_identical(p, before)
  expect_identical(q$fallbacks$sigma$source, "backend_default")
  expect_false(q$fallbacks$sigma$verified)
  expect_false(q$fallbacks$sigma$satisfied)
  expect_identical(q$candidates$sigma$status, "no_supported_fit")
  expect_length(q$beliefs$sigma, 2)
  expect_true(q$backend$compatible)
  text = paste(capture.output(print(summary(q))), collapse = "\n")
  expect_match(text, "Backend default (fallback)", fixed = TRUE)
  expect_match(text, "original rstanarm default has been retained", fixed = TRUE)
  expect_match(text, "does not satisfy all current beliefs", fixed = TRUE)
  expect_false(grepl("None selected", text, fixed = TRUE))
  expect_false(grepl("Finite df search", text, fixed = TRUE))
  details = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
  expect_match(details, "Finite df search", fixed = TRUE)
  plot = plot(q, parameter = "sigma") + ggplot2::theme_bw()
  expect_s3_class(plot, "ggplot")
  expect_match(plot$labels$subtitle, "fallback")
  expect_equal(plot$data$density, dexp(plot$data$value, p$selected$sigma$rate))
  expect_error(as_rstanarm(q), "sigma.*allow_fallback = TRUE")
  expect_error(fit(q), "sigma.*allow_fallback = TRUE")
  expect_error(check_prior(q), "sigma.*allow_fallback = TRUE")
  expect_warning({args = as_rstanarm(q, allow_fallback = TRUE)}, "Using retained fallback")
  expect_identical(args, as_rstanarm(p))
  set.seed(12)
  a = check_prior(p, 100)
  set.seed(12)
  expect_warning({b = check_prior(q, 100, allow_fallback = TRUE)}, "Using retained fallback")
  expect_identical(a$parameter_draws, b$parameter_draws)
  expect_identical(a$expected, b$expected)
  expect_identical(a$simulated, b$simulated)
  expect_true(b$metadata$fallback$approved)
})

test_that("previous elicited priors survive failed additions and recover on removal", {
  p = prior_example() + update_priors("sigma", belief(mean = 30000)) + update_priors("x", belief(mean = 2), belief(sd = 1))
  q = p + update_priors("sigma", belief(mode = 30000))
  expect_identical(q$selected, p$selected)
  expect_identical(q$fallbacks$sigma$source, "elicited")
  expect_identical(q$candidates$sigma$status, "backend_incompatible")
  expect_match(paste(capture.output(print(summary(q))), collapse = "\n"), "Previous elicited (fallback)", fixed = TRUE)
  expect_identical(q$selected$sigma$rate, 1 / 30000)
  restored = update(q, parameter = "sigma", remove = 2)
  expect_null(restored$fallbacks$sigma)
  expect_identical(restored$selected, p$selected)
  expect_no_warning(as_rstanarm(restored))
  replaced = update(q, parameter = "sigma", remove = 1:2, belief(mean = 5000))
  expect_null(replaced$fallbacks$sigma)
  expect_identical(replaced$selected$sigma$rate, 1 / 5000)
  cleared = update(q, parameter = "sigma", remove = 1:2)
  expect_identical(cleared$selected$sigma, p$defaults$sigma)
  expect_identical(cleared$last_valid_source[["sigma"]], "backend_default")
})

test_that("coefficient failures preserve priors and mathematical diagnostics", {
  p = update(prior_example(), parameter = "x", belief(mean = 2), belief(sd = 1))
  q = p + update_priors("x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
  expect_identical(q$selected, p$selected)
  expect_identical(q$candidates$x$status, "contradictory")
  expect_length(q$beliefs$x, 4)
  expect_equal(q$candidates$x$conflicts, c(3, 4))
  expect_s3_class(plot(q), "ggplot")
  expect_s3_class(plot(q, parameter = "x"), "ggplot")
  expect_error(as_rstanarm(q), "contradictory")
  expect_warning({args = as_rstanarm(q, allow_fallback = TRUE)}, "contradictory")
  expect_identical(args, as_rstanarm(p))
  restored = update(q, parameter = "x", remove = 3:4)
  expect_identical(restored$selected, p$selected)
  expect_null(restored$fallbacks$x)
})

test_that("backend-compatible retained collections recover across update orders", {
  base = prior_example()
  for (order in list(c("x", "gb"), c("gb", "x"))) {
    p = fallback_student_beliefs(base, order[1])
    expect_identical(p$selected, base$selected)
    expect_identical(p$candidates[[order[1]]]$status, "backend_incompatible")
    ranks = p$candidates[[order[1]]]$alternatives
    q = fallback_student_beliefs(p, order[2])
    expect_length(q$fallbacks, 0)
    expect_equal(q$candidates[[order[1]]]$alternatives, ranks)
    expect_true(all(vapply(q$selected[c("x", "gb")], function(p) p$distribution == "student_t", logical(1))))
    # Resetting one slope must not silently discard the other slope's valid t.
    cleared = update(q, parameter = order[2], remove = 1:3)
    expect_identical(cleared$selected, q$selected)
    expect_true(cleared$backend$compatible)
    expect_identical(cleared$candidates[[order[2]]]$status, "backend_incompatible")
    expect_error(as_rstanarm(cleared), "allow_fallback")
    expect_warning(as_rstanarm(cleared, allow_fallback = TRUE), "fallback")
    restored = update(cleared, parameter = order[1], remove = 1:3)
    expect_identical(restored$selected, base$selected)
    expect_length(restored$fallbacks, 0)
  }
})

test_that("multiple failures and numerical failures preserve independent successes", {
  p = update(prior_example(), parameter = "x", belief(mean = 2))
  q = p + update_priors("sigma", belief(mode = 5)) +
    update_priors("x", belief(mean = -100))
  expect_identical(q$selected, p$selected)
  expect_setequal(prior_fallback_names(q), c("sigma", "x"))
  expect_error(as_rstanarm(q), "x.*sigma|sigma.*x")
  expect_warning(as_rstanarm(q, allow_fallback = TRUE), "fallback")
  local_mocked_bindings(prior_fit_candidates = function(...) stop("forced solver failure"))
  failed = p + update_priors("x", belief(sd = 3))
  expect_identical(failed$selected, p$selected)
  expect_identical(failed$candidates$x$status, "no_supported_fit")
  expect_match(failed$fallbacks$x$reason, "mathematical incompatibility was not established", fixed = TRUE)
})

test_that("failed outcome reconciliation retains a coherent previous joint collection", {
  p = joint_example() + update_outcome(list(x = 6, z = 1), belief(mean = 12))
  expect_warning({q = p + update_outcome(list(x = 6, z = 1), belief(mean = 100))}, "incompatible")
  expect_identical(q$selected, p$selected)
  expect_true(q$joint$retained)
  expect_identical(q$joint$status, "incompatible")
  expect_length(q$outcomes, 2)
  expect_error(as_rstanarm(q), "outcome beliefs.*incompatible")
  expect_warning({args = as_rstanarm(q, allow_fallback = TRUE)}, "outcome beliefs")
  expect_identical(args, as_rstanarm(p))
  expect_s3_class(plot(q), "ggplot")
  text = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
  expect_match(text, "Available priors (retained", fixed = TRUE)
  expect_match(text, "No verified reconciliation", fixed = TRUE)
  # Sigma remains independent, even when an earlier outcome update failed.
  expect_warning({revised = q + update_priors("sigma", belief(mean = 4))}, "incompatible")
  expect_identical(revised$selected$sigma, list(distribution = "exponential", rate = .25))
  expect_identical(revised$selected[colnames(p$model$X)], p$selected[colnames(p$model$X)])
})

test_that("numerical outcome failures remain unverified and preserve priors", {
  p = joint_example()
  local_mocked_bindings(prior_joint_solve = function(...) stop("forced numerical failure"))
  expect_warning({q = p + update_outcome(list(x = 6, z = 1), belief(mean = 100))}, "No mathematical incompatibility")
  expect_identical(q$selected, p$selected)
  expect_identical(q$joint$status, "unresolved")
  expect_true(q$joint$retained)
  expect_null(q$joint$compatibility)
  expect_error(check_prior(q), "unresolved")
  expect_warning(check_prior(q, 10, allow_fallback = TRUE), "outcome beliefs")
})

test_that("fallback opt-in never bypasses stale data, tampering, or backend checks", {
  q = prior_example() + update_priors("sigma", belief(mode = 10))
  for (allow in list(NA, 1, "yes", logical())) expect_error(as_rstanarm(q, allow_fallback = allow), "allow_fallback must")
  altered = q
  altered$selected$sigma$rate = altered$selected$sigma$rate * 2
  expect_error(as_rstanarm(altered, allow_fallback = TRUE), "Retained prior")
  altered = q
  altered$beliefs$sigma[[1]]$target = 5
  expect_error(check_prior(altered, allow_fallback = TRUE), "metadata changed")
  altered = q
  altered$model$data$y[1] = altered$model$data$y[1] + 1
  expect_error(as_rstanarm(altered, allow_fallback = TRUE), "changed")
  altered = q
  altered$model$reference[1] = altered$model$reference[1] + 1
  expect_error(check_prior(altered, allow_fallback = TRUE), "parameterization")
  altered = q
  altered$backend$compatible = FALSE
  expect_error(check_prior(altered, allow_fallback = TRUE), "incompatible with the backend")
})

test_that("fit consumes opt-in and fully successful priors behave as before", {
  p = prior_example()
  expect_identical(as_rstanarm(p, allow_fallback = TRUE), as_rstanarm(p))
  set.seed(20)
  a = check_prior(p, 10)
  set.seed(20)
  expect_identical(a, check_prior(p, 10, allow_fallback = TRUE))
  q = p + update_priors("sigma", belief(mode = 10))
  captured = new.env()
  local_mocked_bindings(prior_rstanarm_fit = function(arguments) { captured$arguments = arguments; "sentinel" })
  expect_error(fit(q, chains = 1), "allow_fallback")
  expect_warning({result = fit(q, chains = 1, allow_fallback = TRUE)}, "fallback")
  expect_identical(result, "sentinel")
  expect_equal(captured$arguments$chains, 1)
  expect_false("allow_fallback" %in% names(captured$arguments))
  expect_identical(captured$arguments$prior_aux, as_rstanarm(p)$prior_aux)
})
