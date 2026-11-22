prior_removal_example = function(backend) {
  p = prior_example()
  if (backend == "rstanarm") return(p)
  skip_if_not_installed("brms")
  identify_priors(p$model$formula, p$model$data, package = backend)
}

test_that("composable removals preserve positions, order, IDs and unrelated priors", {
  for (backend in c("rstanarm", "brms")) {
    p = prior_removal_example(backend) + update_priors("x", belief(mode = 2), belief(sd = 1),
      belief(value < 3, probability = pnorm(1)), belief(mean = 2)) + update_priors("gb", belief(mean = -1))
    original = p
    one = p + update_priors("x", remove = 2)
    expect_identical(one, prior_apply_parameter(p, parameter = "x", remove = 2))
    expect_identical(one$beliefs$x, p$beliefs$x[-2])
    multiple = p + update_priors("x", remove = c(4, 2))
    expect_identical(multiple, prior_apply_parameter(p, parameter = "x", remove = c(4, 2)))
    expect_identical(multiple$beliefs$x, p$beliefs$x[c(1, 3)])
    expect_identical(multiple$next_id, p$next_id)
    expect_identical(multiple$selected$gb, p$selected$gb)
    expect_identical(multiple$beliefs$gb, p$beliefs$gb)
    expect_identical(multiple$selected$sigma, p$selected$sigma)
    all = p + update_priors("x", remove = 1:4)
    expect_identical(all, prior_apply_parameter(p, parameter = "x", remove = 1:4))
    expect_length(all$beliefs$x, 0)
    expect_identical(all$selected$x, p$defaults$x)
    expect_null(all$fallbacks$x)
    expect_identical(all$selected$gb, p$selected$gb)
    expect_identical(p, original)
  }
})

test_that("replacements remove simultaneously before appending in supplied order", {
  for (backend in c("rstanarm", "brms")) {
    p = prior_removal_example(backend) + update_priors("x", belief(mode = 2), belief(sd = 1),
      belief(mean = 2), belief(value < 3, probability = pnorm(1)))
    spec = update_priors("x", belief(sd = .5), belief(value < 3, probability = pnorm(2)), remove = c(4, 2))
    saved = spec
    q = p + spec
    expect_identical(q, prior_apply_parameter(p, parameter = "x", belief(sd = .5),
      belief(value < 3, probability = pnorm(2)), remove = c(4, 2)))
    expect_identical(q$beliefs$x[1:2], p$beliefs$x[c(1, 3)])
    expect_equal(vapply(q$beliefs$x, `[[`, numeric(1), "id"), c(1, 3, 5, 6))
    expect_equal(vapply(q$beliefs$x, `[[`, character(1), "type"), c("mode", "mean", "sd", "probability"))
    expect_equal(q$selected$x$sd, .5)
    expect_identical(spec, saved)
    expect_identical(q + update_priors("x", remove = c(1, 2)),
      prior_apply_parameter(q, parameter = "x", remove = c(1, 2)))
  }
})

test_that("parameter specifications and their application validate removal indices without partial updates", {
  invalid = list(NULL, NA_real_, NaN, Inf, -Inf, 1 + 1i, "1", TRUE, list(1), 0, -1, .5, c(1, 1))
  for (backend in c("rstanarm", "brms")) {
    p = prior_removal_example(backend) + update_priors("x", belief(mean = 2))
    original = p
    for (remove in invalid) {
      expect_error(update_priors("x", remove = remove), "Invalid removal numbers")
      expect_error(prior_apply_parameter(p, parameter = "x", remove = remove), "Invalid removal numbers")
      expect_error(update_priors("x", belief(mean = 3), remove = remove), "Invalid removal numbers")
    }
    spec = update_priors("x", belief(mean = 3), remove = 2)
    expect_error(p + spec, "Invalid removal numbers.*existing belief list")
    expect_error(prior_apply_parameter(p, parameter = "x", belief(mean = 3), remove = 2), "Invalid removal numbers.*existing belief list")
    expect_error(p + update_priors("sigma", remove = 1), "Invalid removal numbers")
    expect_error(p + update_priors("missing", remove = 1), "Unknown parameter")
    expect_error(update_priors("x", remove = belief(mean = 1)), "Invalid removal numbers")
    expect_error(update_priors("x", 1, remove = 1), "belief objects")
    expect_error(p + update_priors("x", belief(mean = 1, scale = "probability"), remove = 1), "scale|intercept|probability")
    expect_identical(p, original)
    expect_identical(p + update_priors("x", belief(sd = 1), remove = integer()),
      prior_apply_parameter(p, parameter = "x", belief(sd = 1)))

  }
})

test_that("removing failed beliefs clears fallback safeguards through the shared parameter update", {
  for (backend in c("rstanarm", "brms")) {
    p = prior_removal_example(backend) + update_priors("x", belief(mean = 2), belief(sd = 1))
    bad = p + update_priors("x", belief(value < 0, probability = .95), belief(value > 1, probability = .95))
    expect_identical(bad$candidates$x$status, "contradictory")
    expect_error(check_prior(bad), "allow_fallback")
    restored = bad + update_priors("x", remove = c(3, 4))
    expect_identical(restored, prior_apply_parameter(bad, parameter = "x", remove = c(3, 4)))
    expect_null(restored$fallbacks$x)
    expect_identical(restored$selected, p$selected)
    expect_s3_class(check_prior(restored, 5), "flexbayes_prior_predictive")
    export = if (backend == "brms") as_brms else as_rstanarm
    expect_type(export(restored), "list")
    unsupported = p + update_priors("sigma", belief(mean = 2), belief(mean = 3))
    expect_identical(unsupported$candidates$sigma$status, "no_supported_fit")
    resolved = unsupported + update_priors("sigma", belief(mean = 4), remove = c(1, 2))
    expect_null(resolved$fallbacks$sigma)
    expect_identical(resolved$selected$sigma$rate, .25)
    expect_identical(resolved, prior_apply_parameter(unsupported, parameter = "sigma", belief(mean = 4), remove = c(1, 2)))
  }
})

test_that("removals preserve whole-model compatibility and outcome reconciliation", {
  for (backend in c("rstanarm", "brms")) {
    p = prior_removal_example(backend)
    t_beliefs = list(belief(mean = 0), belief(sd = sqrt(7 / 5)), belief(value < 1, probability = pt(1, 7)))
    p = (p + do.call(update_priors, c(list(parameter = "x"), t_beliefs)))
    p = (p + do.call(update_priors, c(list(parameter = "gb"), t_beliefs)))
    reset = p + update_priors("x", remove = 1:3)
    expect_identical(reset, prior_apply_parameter(p, parameter = "x", remove = 1:3))
    expect_identical(reset$selected$gb, p$selected$gb)
    if (backend == "rstanarm") {
      # Reset is deferred rather than changing another slope's valid t prior.
      expect_identical(reset$fallbacks$x$status, "backend_incompatible")
      expect_error(as_rstanarm(reset), "allow_fallback")
    } else expect_identical(reset$selected$x, p$defaults$x)
    all = reset + update_priors("gb", remove = 1:3)
    expect_identical(all$selected$x, p$defaults$x)
    expect_identical(all$selected$gb, p$defaults$gb)
    expect_length(all$fallbacks, 0)
    d = data.frame(x = 1:6, z = rep(0:1, 3), y = c(2, 4, 3, 7, 6, 10))
    joint = identify_priors(y ~ x + z, d, package = backend) + update_priors("x", belief(mean = 2), belief(sd = 1)) +
      update_priors(outcome = list(x = 3, z = 0), belief(mean = 6))
    q = joint + update_priors("x", belief(sd = 2), remove = 2)
    expect_identical(q, prior_apply_parameter(joint, parameter = "x", belief(sd = 2), remove = 2))
    expect_identical(q$outcomes, joint$outcomes)
    expect_identical(q$selected$sigma, joint$selected$sigma)
  }
})
