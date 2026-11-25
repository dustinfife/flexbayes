joint_numeric_cache = new.env(parent = emptyenv())

joint_housing_base = function(package = "rstanarm", response_scale = 1, sqft_scale = 1) {
  skip_if_not_installed(package)
  d = homes_2024
  d$sale_price = d$sale_price * response_scale
  d$sqft = d$sqft * sqft_scale
  identify_priors(sale_price ~ age + sqft + bedrooms, d, package = package) +
    update_priors("(Intercept)", belief(mean = 417000 * response_scale, tolerance = 1e-6 * response_scale),
      belief(value < 300000 * response_scale, probability = .05)) +
    update_priors("age", belief(mode = -.2, scale = "standardized"),
      belief(value > -.01, probability = .025, scale = "standardized")) +
    update_priors("bedrooms", belief(mode = 17097 * response_scale, tolerance = 1e-6 * response_scale),
      belief(value < 0, probability = .01))
}

joint_housing_update = function(p, probability = TRUE, second = TRUE, response_scale = 1, sqft_scale = 1) {
  records = list(belief(mode = 430000 * response_scale, tolerance = 1000 * response_scale))
  if (probability) records = c(records, list(belief(value < 400000 * response_scale,
    probability = .30, tolerance = .005)))
  p = p + do.call(update_priors, c(list(outcome = list(sqft = 1700 * sqft_scale, age = 44, bedrooms = 3)), records))
  if (!second) return(p)
  p + update_priors(outcome = list(sqft = 2200 * sqft_scale, age = 44, bedrooms = 3),
    belief(mode = 530000 * response_scale, tolerance = 1000 * response_scale))
}

joint_housing_fixture = function(package = "rstanarm") {
  if (exists(package, joint_numeric_cache, inherits = FALSE)) return(get(package, joint_numeric_cache))
  p = joint_housing_base(package)
  q = suppressWarnings(joint_housing_update(p))
  result = list(base = p, result = q)
  assign(package, result, joint_numeric_cache)
  result
}

joint_numeric_moments = function(p, at) {
  X = model.matrix(delete.response(terms(p$model$formula)), as.data.frame(at))
  columns = setdiff(colnames(X), "(Intercept)")
  X[, columns] = sweep(X[, columns, drop = FALSE], 2, colMeans(p$model$X[, columns, drop = FALSE]), "-")
  mean = vapply(p$selected[colnames(X)], `[[`, numeric(1), "mean")
  sd = vapply(p$selected[colnames(X)], `[[`, numeric(1), "sd")
  c(mean = sum(X * mean), sd = sqrt(sum((X * sd)^2)))
}

joint_numeric_jacobian = function(fun, x, step = 1e-5) {
  vapply(seq_along(x), function(i) {
    above = below = x
    above[i] = above[i] + step
    below[i] = below[i] - step
    (fun(above) - fun(below)) / (2 * step)
  }, numeric(length(fun(x))))
}

test_that("the exact housing example verifies all three tiers for both packages", {
  for (package in c("rstanarm", "brms")) {
    fixture = joint_housing_fixture(package)
    p = fixture$base
    q = fixture$result
    expect_identical(q$joint$status, "reconciled")
    expect_length(q$fallbacks, 0)
    expect_false(isTRUE(q$joint$retained))
    expect_identical(q$beliefs, p$beliefs)
    expect_identical(q$defaults, p$defaults)
    expect_identical(q$selected$sigma, p$selected$sigma)
    expect_match(q$joint$stages$KL$method, "Newton barrier", fixed = TRUE)
    for (stage in q$joint$stages) {
      expect_lt(stage$violation, 1e-8)
      expect_lt(stage$stationarity, 1e-4)
      expect_true(is.finite(stage$objective))
    }
    moments = joint_numeric_moments(q, list(age = 44, sqft = 1700, bedrooms = 3))
    other = joint_numeric_moments(q, list(age = 44, sqft = 2200, bedrooms = 3))
    expect_lte(abs(moments["mean"] - 430000), 1000)
    expect_lte(abs(pnorm(400000, moments[1], moments[2]) - .30), .005)
    expect_lte(abs(other["mean"] - 530000), 1000)
    expect_equal(q$selected$sqft$mean, unname((other[1] - moments[1]) / 500), tolerance = 1e-10)
    expect_gte(q$selected$sqft$mean, 196)
    expect_lte(q$selected$sqft$mean, 204)
    expect_equal(round((400000 - 430000) / qnorm(.30)), 57208)
    expect_equal(pnorm(400000, moments[1], moments[2]), q$joint$outcome$achieved[2], tolerance = 1e-12)
    parameters = prior_joint_parameters(q$joint$theta, q$joint$problem)
    expect_equal(vapply(q$selected[1:4], `[[`, numeric(1), "mean"), parameters$mean)
    expect_equal(vapply(q$selected[1:4], `[[`, numeric(1), "sd"), parameters$sd)
    reference = q$joint$problem
    kl = sum(((parameters$mean - reference$mean) / reference$sd)^2 +
      (parameters$sd / reference$sd)^2 - 1 - 2 * log(parameters$sd / reference$sd)) / 2
    expect_equal(q$joint$kl, kl, tolerance = 1e-12)
    expect_lte(kl, prior_joint_kl(head(q$joint$stages$coefficient$par, 8), reference)$value + 1e-8)
    limit = q$joint$coefficient_limit * (1 + 1e-6)
    expect_lte(max(q$joint$coefficient$normalized_violation), limit + 1e-10)

    # Independently solve the active KKT equations, including nonnegative duals.
    system = prior_joint_system(q$joint$theta, reference, 3, limit)
    active = which(system$g >= -1e-7)
    multipliers = qr.solve(t(system$jacobian[active, , drop = FALSE]), -system$gradient, tol = 1e-12)
    expect_gte(min(multipliers), -1e-8)
    expect_lt(max(abs(system$gradient + crossprod(system$jacobian[active, , drop = FALSE], multipliers))), 1e-4)
    expect_lt(max(abs(multipliers * system$g[active])), 1e-4)

    set.seed(117)
    check = check_prior(q, 100, newdata = data.frame(age = 44, sqft = 1700, bedrooms = 3))
    expect_identical(check$metadata$priors, q$selected)
    expect_equal(check$expected, check$coefficient_draws %*% t(check$X))
    arguments = if (package == "brms") as_brms(q) else as_rstanarm(q)
    expect_type(arguments, "list")
    if (package == "rstanarm") {
      expect_equal(unname(arguments$prior$location), unname(parameters$mean[-1]))
      expect_equal(unname(arguments$prior$scale), unname(parameters$sd[-1]))
    } else expect_match(arguments$prior$prior[arguments$prior$coef == "sqft"], "normal(", fixed = TRUE)
  }
})

test_that("location-only, location/probability, and multiple-location variants remain verified", {
  for (package in c("rstanarm", "brms")) {
    p = joint_housing_fixture(package)$base
    for (probability in c(FALSE, TRUE)) for (second in c(FALSE, TRUE)) {
      q = suppressWarnings(joint_housing_update(p, probability, second))
      expect_true(q$joint$status %in% c("compatible", "reconciled"))
      expect_true(all(q$joint$outcome$satisfied))
      expect_length(q$fallbacks, 0)
    }
    q = joint_housing_fixture(package)$result
    at = list(age = 55, sqft = 2000, bedrooms = 4)
    target = unname(joint_numeric_moments(q, at)[1])
    q = suppressWarnings(q + update_priors(outcome = at, belief(mean = target, tolerance = 1000)))
    expect_true(all(q$joint$outcome$satisfied))
    expect_true(q$joint$status %in% c("compatible", "reconciled"))
    bad = suppressWarnings(q + update_priors(outcome = list(age = 44, sqft = 1700, bedrooms = 3), belief(mean = 1)))
    expect_identical(bad$joint$status, "incompatible")
    expect_identical(bad$selected, q$selected)
    expect_error(check_prior(bad), "allow_fallback")
  }
})

test_that("the recovered solution is deterministic and invariant to equivalent units", {
  fixture = joint_housing_fixture()
  set.seed(95)
  before = .Random.seed
  repeated = suppressWarnings(joint_housing_update(fixture$base))
  expect_identical(.Random.seed, before)
  expect_equal(repeated$selected, fixture$result$selected, tolerance = 1e-10)
  p = joint_housing_base(response_scale = .001, sqft_scale = .001)
  q = suppressWarnings(joint_housing_update(p, response_scale = .001, sqft_scale = .001))
  expect_true(all(q$joint$outcome$satisfied))
  expect_true(q$joint$status %in% c("compatible", "reconciled"))
  conversions = c(.001, .001, 1, .001)
  base = fixture$result
  means = vapply(q$selected[1:4], `[[`, numeric(1), "mean") / conversions
  sd = vapply(q$selected[1:4], `[[`, numeric(1), "sd") / conversions
  original_sd = vapply(fixture$base$selected[1:4], `[[`, numeric(1), "sd")
  expect_equal(means / original_sd, vapply(base$selected[1:4], `[[`, numeric(1), "mean") / original_sd, tolerance = 2e-6)
  expect_equal(sd / original_sd, vapply(base$selected[1:4], `[[`, numeric(1), "sd") / original_sd, tolerance = 2e-6)
  expect_equal(q$joint$coefficient_limit, base$joint$coefficient_limit, tolerance = 2e-6)
})

test_that("a convergence code never overrides failed verification or misleading diagnostics", {
  p = joint_example(y ~ 1) + update_priors("(Intercept)", belief(mean = 0), belief(sd = 1))
  local_mocked_bindings(prior_joint_solve = function(start, ...) list(par = start, convergence = 0L,
    iterations = 1L, objective = 0, violation = 0, stationarity = 1))
  expect_warning({q = p + update_priors(outcome = list(), belief(mean = 10))}, "could not verify")
  expect_identical(q$joint$status, "unresolved")
  expect_identical(q$selected, p$selected)
  expect_identical(q$beliefs, p$beliefs)
  expect_error(check_prior(q), "allow_fallback")
  advanced = paste(capture.output(print(summary(q, diagnostics = TRUE))), collapse = "\n")
  expect_match(advanced, "Original coefficient prior:", fixed = TRUE)
  expect_match(advanced, "Retained fallback prior:", fixed = TRUE)
  expect_match(advanced, "Attempted numerical prior (not selected):", fixed = TRUE)
  expect_match(advanced, prior_format(list(distribution = "normal", mean = 10, sd = 1)), fixed = TRUE)
  expect_match(advanced, "Attempted change in mean: 10", fixed = TRUE)
  expect_match(advanced, "KL(attempted || baseline)", fixed = TRUE)
  ordinary = paste(capture.output(print(summary(q))), collapse = "\n")
  expect_false(grepl("Attempted numerical", ordinary, fixed = TRUE))
  expect_false(grepl("stationarity", ordinary, fixed = TRUE))
  expect_warning(as_rstanarm(q, allow_fallback = TRUE), "Using retained fallback")
})

test_that("equivalent face constraints have correct gradients and Hessians", {
  q = joint_housing_fixture()$result
  problem = q$joint$problem
  face = list(problem = problem, sd = prior_joint_parameters(head(q$joint$stages$coefficient$par, 8), problem)$sd,
    limit = q$joint$coefficient_limit * (1 + 1e-6))
  x = c(head(q$joint$theta, 4), vapply(q$selected[1:4], `[[`, numeric(1), "sd") / face$sd)
  system = prior_joint_face_system(x, face)
  near = c(rep(1e-8, 4), (1 + 1e-8) * problem$sd / face$sd)
  value = prior_joint_face_system(near, face)$value
  reference = prior_joint_kl(c(rep(1e-8, 4), rep(log1p(1e-8), 4)), problem)$value
  expect_gte(value, 0)
  expect_lt(abs(value - reference), 1e-21)
  expect_lt(max(abs(joint_numeric_jacobian(function(z) prior_joint_face_system(z, face)$value, x) - system$gradient)), 1e-6)
  expect_lt(max(abs(joint_numeric_jacobian(function(z) prior_joint_face_system(z, face)$g, x) - system$jacobian)), 1e-6)
  expect_equal(joint_numeric_jacobian(function(z) prior_joint_face_system(z, face)$gradient, x), system$hessian, tolerance = 1e-5, ignore_attr = TRUE)
  derivatives = array(joint_numeric_jacobian(function(z) as.vector(prior_joint_face_system(z, face)$jacobian), x),
    c(length(system$g), 8, 8))
  for (i in seq_along(system$g)) expect_lt(max(abs(derivatives[i, , ] - system$curvature[[i]])), 1e-6)
  broad = face
  broad$limit = 3
  for (i in seq_along(broad$problem$outcome)) broad$problem$outcome[[i]]$tolerance = broad$problem$outcome[[i]]$tolerance * 20
  for (phase in c(FALSE, TRUE)) {
    point = if (phase) c(x, 1) else x
    barrier = prior_joint_face_barrier(point, broad, .01, phase)
    expect_true(is.finite(barrier$value))
    expect_equal(as.vector(joint_numeric_jacobian(function(z) prior_joint_face_barrier(z, broad, .01, phase)$value, point)),
      unname(barrier$gradient), tolerance = 1e-5)
    expect_equal(joint_numeric_jacobian(function(z) prior_joint_face_barrier(z, broad, .01, phase)$gradient, point),
      barrier$hessian, tolerance = 1e-5, ignore_attr = TRUE)
  }
})

test_that("the housing minimax value agrees with an independent variance lower bound", {
  fixture = joint_housing_fixture()
  p = fixture$base
  q = fixture$result
  problem = q$joint$problem
  row = c(1, 44 - mean(p$model$data$age), 1700 - mean(p$model$data$sqft), 3 - mean(p$model$data$bedrooms))
  first = 430000 + c(-1, 1) * 1000 * (1 - problem$interior)
  second = 530000 - 1000 * (1 - problem$interior)
  probability = .30 + .005 * (1 - problem$interior)
  baseline_sd = vapply(p$selected[1:4], `[[`, numeric(1), "sd")
  baseline_mean = vapply(p$selected[1:4], `[[`, numeric(1), "mean")
  lower_bound = function(t) {
    # For this design, increasing age and decreasing bedrooms minimizes both
    # the required intercept and their own minimum admissible variances.
    age = baseline_mean[2] + p$beliefs$age[[1]]$tolerance + t * baseline_sd[2]
    bedrooms = baseline_mean[4] - p$beliefs$bedrooms[[1]]$tolerance - t * baseline_sd[4]
    vapply(first, function(location) {
      sqft = (second - location) / 500
      intercept = location - sum(row[-1] * c(age, sqft, bedrooms))
      sd = c((intercept - 300000) / (-qnorm(.05 - 1e-6) + t),
        (p$beliefs$age[[2]]$threshold - age) / (-qnorm(.025 - 1e-6) + t), 0,
        bedrooms / (-qnorm(.01 - 1e-6) + t))
      sum((row * sd)^2) - ((location - 400000) / qnorm(probability))^2
    }, numeric(1))
  }
  # The variance deficit is concave in the first location, so the smallest
  # possible deficit occurs at an endpoint, not an unexamined interior point.
  t = q$joint$coefficient_limit
  expect_lt(((1 + row[3] / 500) / (-qnorm(.05 - 1e-6) + t))^2 - 1 / qnorm(probability)^2, 0)
  root = uniroot(function(t) min(lower_bound(t)), c(1.5, 1.7), tol = 1e-13)$root
  expect_gte(t, root - 1e-10)
  expect_lt((t - root) / root, 1e-8)
  expect_true(all(vapply(q$selected[1:4], `[[`, numeric(1), "sd") > 0))
  expect_gt(q$selected$sqft$sd, 0)
})

test_that("minimax epigraph and near-degenerate face derivatives remain faithful", {
  q = joint_housing_fixture()$result
  problem = q$joint$problem
  for (i in seq_along(problem$outcome)) problem$outcome[[i]]$tolerance = problem$outcome[[i]]$tolerance * 20
  face = list(problem = problem, sd = prior_joint_parameters(q$joint$theta, problem)$sd,
    limit = 3, minimax = TRUE, coefficient_scale = 2)
  x = c(head(q$joint$theta, 4), rep(1, 4), 3)
  system = prior_joint_minimax_system(x, face)
  expect_lt(max(abs(joint_numeric_jacobian(function(z) prior_joint_minimax_system(z, face)$g, x) - system$jacobian)), 1e-6)
  derivative = array(joint_numeric_jacobian(function(z) as.vector(prior_joint_minimax_system(z, face)$jacobian), x),
    c(length(system$g), 9, 9))
  for (i in seq_along(system$g)) expect_lt(max(abs(derivative[i, , ] - system$curvature[[i]])), 1e-6)
  for (phase in c(FALSE, TRUE)) {
    point = if (phase) c(x, 1) else x
    barrier = prior_joint_face_barrier(point, face, .01, phase)
    expect_equal(as.vector(joint_numeric_jacobian(function(z) prior_joint_face_barrier(z, face, .01, phase)$value, point)),
      unname(barrier$gradient), tolerance = 1e-5)
    expect_equal(joint_numeric_jacobian(function(z) prior_joint_face_barrier(z, face, .01, phase)$gradient, point),
      barrier$hessian, tolerance = 1e-5, ignore_attr = TRUE)
  }
  p = joint_example(y ~ 1) + update_priors("(Intercept)", belief(mean = 0), belief(sd = 1))
  p = suppressWarnings(p + update_priors(outcome = list(), belief(mean = 10)))
  face = list(problem = p$joint$problem, sd = 1, limit = 11)
  for (spread in c(1e-6, 1, 1000)) {
    system = prior_joint_face_system(c(10, spread), face)
    expect_true(all(is.finite(c(system$value, system$gradient, system$jacobian))))
    expect_true(all(vapply(system$curvature, function(H) all(is.finite(H)), logical(1))))
  }
  expect_identical(prior_joint_face_barrier(c(10, 0), face, .01)$value, Inf)
  for (operator in c("<", ">")) for (probability in c(0, 1)) {
    threshold = if ((operator == "<") == (probability == 0)) -1e6 else 1e6
    record = if (operator == "<") belief(value < threshold, probability = probability) else
      belief(value > threshold, probability = probability)
    record$configuration = 1L
    face$problem$outcome = list(record)
    original = prior_joint_system(c(10, 0), face$problem, 3, 11)
    equivalent = prior_joint_face_system(c(10, 1), face)
    expect_true(all(is.finite(c(equivalent$g, equivalent$jacobian))))
    expect_identical(all(equivalent$g <= 0), all(original$g <= 0))
  }
})
