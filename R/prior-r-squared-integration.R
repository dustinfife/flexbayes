# Numerical limits are computation limits, never a relaxation of belief bands.
prior_r_squared_policy = function() list(orders = c(16L, 32L, 64L, 128L),
  maximum_dimension = 3L, maximum_points = 4500000L, tail = 1e-10)
prior_is_r_squared = function(record) identical(record$scale, "r_squared")
prior_has_r_squared = function(object) any(vapply(object$beliefs$sigma, prior_is_r_squared, logical(1)))

prior_predictor_covariance = function(X) {
  slopes = setdiff(colnames(X), "(Intercept)")
  if (!length(slopes)) return(matrix(numeric(), 0, 0, dimnames = list(character(), character())))
  stats::cov(X[, slopes, drop = FALSE])
}

prior_r_squared_belief = function(record, object, parameter) {
  if (!prior_has_sigma(object) || object$model$family$link != "identity" || parameter != "sigma") stop("r_squared beliefs require sigma in a Gaussian identity model.")
  if (record$type != "probability") stop("r_squared currently supports probability/quantile beliefs only, not mean, mode or SD beliefs.")
  if (record$threshold < 0 || record$threshold > 1) stop("R-squared thresholds must be between 0 and 1.")
  record
}

prior_r_squared_context = function(object) {
  if (!prior_has_sigma(object) || object$model$family$link != "identity") stop("R-squared elicitation requires a Gaussian identity model.")
  prior_validate_standardization_data(object)
  covariance = prior_predictor_covariance(object$model$X)
  if (is.null(object$model$predictor_covariance) ||
      !isTRUE(all.equal(covariance, object$model$predictor_covariance))) stop("Recorded predictor covariance changed; identify priors again.")
  if (any(!is.finite(covariance))) stop("Finite predictor covariance required for R-squared elicitation.")
  slopes = colnames(covariance)
  priors = object$selected[slopes]
  for (prior in priors) {
    prior_validate_distribution(prior)
    if (!prior$distribution %in% c("normal", "student_t")) stop("Unsupported coefficient prior for R-squared integration.")
  }
  active = which(diag(covariance) > 0)
  covariance_active = covariance[active, active, drop = FALSE]
  priors_active = priors[active]
  scale = vapply(priors_active, function(p) if (p$distribution == "normal") p$sd else p$scale, numeric(1))
  location = vapply(priors_active, function(p) if (p$distribution == "normal") p$mean else p$location, numeric(1))
  if (any(!is.finite(location / scale))) stop("Coefficient standardization exceeds the R-squared integration range.")
  B = covariance_active * outer(scale, scale)
  if (length(B) && (any(!is.finite(B)) || !any(B > 0))) stop("Explained-variance scaling exceeds numerical range.")
  if (length(B) && min(eigen(B, symmetric = TRUE, only.values = TRUE)$values) < -1e-12 * max(abs(B))) stop("Predictor covariance is not positive semidefinite at numerical precision.")
  list(covariance = covariance, coefficients = priors, active = active, B = B,
    center = location / scale, priors = priors_active,
    observed_variance = stats::var(stats::model.response(object$model$frame)))
}

# Gauss-Legendre nodes and weights use only base R's symmetric eigensolver.
prior_r_squared_rule = function(n, lower, upper) {
  j = seq_len(n - 1L)
  off = j / sqrt(4 * j^2 - 1)
  J = matrix(0, n, n)
  J[cbind(j, j + 1L)] = off
  J[cbind(j + 1L, j)] = off
  e = eigen(J, symmetric = TRUE)
  list(x = (upper + lower) / 2 + (upper - lower) / 2 * e$values,
    w = (upper - lower) * e$vectors[1, ]^2)
}

prior_r_squared_density = function(z, prior, center) {
  if (prior$distribution == "normal") stats::dnorm(z - center) else stats::dt(z - center, prior$df)
}

# Integrate in spherical coordinates about zero explained effect, removing the
# norm cusp. Tail truncation uses marginal union bounds, not outcome variance.
prior_r_squared_grid = function(context, order) {
  policy = prior_r_squared_policy()
  d = length(context$priors)
  if (!d) return(list(v = 0, w = 1, tail = 0, mass_error = 0, order = order, zero = TRUE))
  fixed = prior_r_squared_fixed_grid(context)
  if (!is.null(fixed)) return(fixed)
  if (d > policy$maximum_dimension) stop("Deterministic R-squared quadrature currently verifies at most three varying model-matrix coefficient columns; higher dimensions remain numerically unresolved.")
  tail = policy$tail
  bounds = vapply(context$priors, function(p) if (p$distribution == "normal")
    stats::qnorm(tail / (2 * d), lower.tail = FALSE) else stats::qt(tail / (2 * d), p$df, lower.tail = FALSE), numeric(1))
  radius = sqrt(sum(context$center^2)) + sqrt(sum(bounds^2))
  if (!is.finite(radius)) stop("Coefficient quadrature radius exceeds finite numerical range.")
  radial_order = order * if (any(vapply(context$priors, function(p) p$distribution == "student_t", logical(1)))) 2L else 1L
  r = prior_r_squared_rule(radial_order, 0, radius)
  if (d == 1L) {
    z = c(r$x, -r$x)
    w = rep(r$w, 2) * prior_r_squared_density(z, context$priors[[1]], context$center[1])
    return(list(v = context$B[1, 1] * z^2, w = w, tail = tail,
      mass_error = abs(sum(w) - 1), order = order))
  }
  phi = 2 * pi * (seq_len(2L * order) - .5) / (2L * order)
  if (d == 2L) {
    direction = cbind(cos(phi), sin(phi))
    angular_weight = rep(pi / order, length(phi))
  } else {
    angle = prior_r_squared_rule(order, -1, 1)
    indices = expand.grid(phi = seq_along(phi), angle = seq_along(angle$x))
    t = angle$x[indices$angle]
    direction = cbind(sqrt(1 - t^2) * cos(phi[indices$phi]),
      sqrt(1 - t^2) * sin(phi[indices$phi]), t)
    angular_weight = angle$w[indices$angle] * pi / order
  }
  if (length(r$x) * nrow(direction) > policy$maximum_points) stop("R-squared quadrature point budget exceeded; no probability has been verified.")
  quadratic = pmax(0, rowSums((direction %*% context$B) * direction))
  indices = expand.grid(radius = seq_along(r$x), direction = seq_len(nrow(direction)))
  radial = r$x[indices$radius]
  w = r$w[indices$radius] * angular_weight[indices$direction] * radial^(d - 1L)
  for (j in seq_len(d)) w = w * prior_r_squared_density(radial * direction[indices$direction, j], context$priors[[j]], context$center[j])
  list(v = radial^2 * quadratic[indices$direction], w = w, tail = tail,
    mass_error = abs(sum(w) - 1), order = order)
}

prior_r_squared_constant = function(record, zero) {
  if (zero) return(prior_joint_constant_value(record))
  h = record$threshold
  lower = record$operator %in% c("<", "<=")
  if (h == 0) return(as.numeric(!lower))
  if (h == 1) return(as.numeric(lower))
  NULL
}

# Conditional sigma integration is analytic: R2 > h iff sigma < sqrt(V(1-h)/h).
# Fixed common coefficient nodes retain the entire coefficient vector.
prior_r_squared_probability = function(record, prior, grid) {
  constant = prior_r_squared_constant(record, isTRUE(grid$zero))
  if (!is.null(constant)) return(constant)
  if (any(!is.finite(grid$v)) || any(!is.finite(grid$w)) || any(grid$v < 0)) stop("Explained variance or integration weights exceed finite numerical range.")
  if (!is.null(grid$budget)) {
    grid$budget$points = grid$budget$points + length(grid$v)
    if (grid$budget$points > 2e8) stop("R-squared probability-evaluation budget exhausted; numerical verification remains unresolved.")
  }
  threshold = sqrt(grid$v) * sqrt(1 - record$threshold) / sqrt(record$threshold)
  operation = if (record$operator %in% c(">", ">=")) "cdf" else "survival"
  sum(grid$w * prior_distribution_operation(prior, operation, threshold))
}

prior_r_squared_values = function(prior, records, grid) {
  vapply(records, function(record) if (prior_is_r_squared(record))
    prior_r_squared_probability(record, prior, grid) else prior_distribution_target(record, prior), numeric(1))
}

# Success requires the original band to contain the whole estimated error band.
# Successive-order differences estimate quadrature error; the tail bound is an
# explicit upper bound. These estimates are not a formal global error proof.
prior_r_squared_verify = function(prior, records, context) {
  policy = prior_r_squared_policy()
  previous = NULL
  last_difference = Inf
  tolerance = vapply(records, `[[`, numeric(1), "tolerance")
  targets = vapply(records, `[[`, numeric(1), "target")
  r2 = vapply(records, prior_is_r_squared, logical(1))
  for (order in policy$orders) {
    grid = prior_r_squared_cached_grid(context, order)
    values = prior_r_squared_values(prior, records, grid)
    if (!is.null(previous)) {
      difference = abs(values - previous)
      error = ifelse(r2, 4 * difference + grid$tail + grid$mass_error, 0)
      for (i in which(r2)) {
        constant = prior_r_squared_constant(records[[i]], isTRUE(grid$zero))
        if (!is.null(constant)) error[i] = 0
        if (!is.null(grid$bounds) && is.null(constant)) {
          envelope = vapply(grid$bounds, function(v) prior_r_squared_probability(records[[i]], prior,
            list(v = v, w = 1)), numeric(1))
          error[i] = max(abs(envelope - values[i])) + grid$tail
        }
      }
      # Two refinements prevent accepting one accidentally matching coarse pair.
      stable = order >= policy$orders[3] && all(difference[r2] <= last_difference[r2] + 1e-14)
      precise = stable && all(error[r2] <= tolerance[r2] / 8)
      if (precise) return(list(valid = all(abs(values - targets) + error <= tolerance),
        achieved = values, error = error, order = grid$order,
        method = if (!is.null(grid$bounds)) "Coefficient tail-envelope (effectively fixed limit)" else "Deterministic spherical quadrature",
        mass_error = grid$mass_error, tail_bound = grid$tail, precise = TRUE))
      last_difference = difference
    }
    previous = values
  }
  list(valid = FALSE, achieved = values, error = error, order = order,
    method = "Numerically unresolved deterministic quadrature", mass_error = grid$mass_error, tail_bound = grid$tail, precise = FALSE)
}

# A fixed-effect limit is used only with an explicit probability envelope from
# coefficient tail bounds. Coefficient uncertainty is never silently dropped.
prior_r_squared_fixed_grid = function(context) {
  d = length(context$priors)
  if (!d) return(NULL)
  tail = prior_r_squared_policy()$tail
  q = vapply(context$priors, function(p) if (p$distribution == "normal")
    stats::qnorm(tail / (2 * d), lower.tail = FALSE) else stats::qt(tail / (2 * d), p$df, lower.tail = FALSE), numeric(1))
  value = sum(context$B * outer(context$center, context$center))
  radius = sqrt(max(0, max(eigen(context$B, symmetric = TRUE, only.values = TRUE)$values)) * sum(q^2))
  if (!is.finite(value) || value <= 0 || !is.finite(radius) || radius > sqrt(value) * 1e-7) return(NULL)
  list(v = value, w = 1, tail = tail, mass_error = 0, order = 0L,
    bounds = c(max(0, sqrt(value) - radius)^2, (sqrt(value) + radius)^2))
}
