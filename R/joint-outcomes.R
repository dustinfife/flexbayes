# Validate an expected-response specification before joint reconciliation.
prior_outcome_update = function(at, records, target = "outcome") {
  # Validate the configuration without changing any of the supplied answers.
  if (!is.list(at) || length(at) && (is.null(names(at)) || anyNA(names(at)) ||
      any(!nzchar(names(at))) || anyDuplicated(names(at)))) {
    stop(target, " must be a uniquely named list.")
  }
  if (any(!vapply(at, function(x) is.numeric(x) && !is.complex(x) && length(x) == 1L &&
      is.null(dim(x)) && is.finite(x), logical(1)))) {
    stop("Each predictor configuration value in ", target, " must be one finite number.")
  }

  # Outcome beliefs concern expected responses, not individual observations.
  if (!length(records) || any(!vapply(records, inherits, logical(1), "prior_belief"))) stop("Supply belief() objects.")
  if (any(!vapply(records, function(r) r$scale == "original" &&
      r$type %in% c("mean", "mode", "probability"), logical(1)))) {
    stop("Expected-outcome elicitation supports original-scale mean, mode and probability beliefs only.")
  }
  list(at = at, beliefs = unname(records))
}

prior_apply_outcome = function(object, specification) {
  for (i in seq_along(specification$beliefs)) {
    specification$beliefs[[i]]$id = object$next_id
    object$next_id = object$next_id + 1L
  }
  record = if (identical(specification$target, "difference"))
    list(difference = specification$difference, beliefs = specification$beliefs) else
    list(at = specification$at, beliefs = specification$beliefs)
  object$outcomes = c(object$outcomes, list(record))
  prior_select_compatible(object)
}

prior_joint_design = function(object, specifications) {
  if (object$model$family$family != "gaussian" || object$model$family$link != "identity") {
    stop("Joint outcome elicitation currently supports Gaussian identity models only.")
  }
  if (!"(Intercept)" %in% colnames(object$model$X)) stop("Joint outcome elicitation requires an intercept.")

  # Reuse the original model matrix and its centering convention.
  variables = all.vars(stats::delete.response(object$model$terms))
  response = attr(object$model$terms, "response")
  predictors = if (response) object$model$frame[-response] else object$model$frame
  if (any(!vapply(object$model$data[variables], is.numeric, logical(1))) ||
      any(!vapply(predictors, is.numeric, logical(1)))) {
    stop("Joint outcome elicitation currently requires numeric predictors, including their formula transformations.")
  }

  rows = lapply(specifications, function(specification) {
    if (!is.null(specification$difference)) {
      from = prior_joint_configuration(object, specification$difference$from, variables, "difference$from")
      to = prior_joint_configuration(object, specification$difference$to, variables, "difference$to")
      return(to - from)
    }
    prior_joint_configuration(object, specification$at, variables, "outcome")
  })
  do.call(rbind, rows)
}

# Parameters are dimensionless mean displacements and log-SD ratios.
prior_joint_parameters = function(theta, problem) {
  n = length(problem$mean)
  list(mean = problem$mean + problem$sd * theta[seq_len(n)],
    sd = problem$sd * exp(theta[n + seq_len(n)]))
}

prior_joint_values = function(theta, problem, coefficient = FALSE) {
  parameters = prior_joint_parameters(theta, problem)
  records = if (coefficient) problem$coefficient else problem$outcome
  p = length(parameters$mean)
  value = numeric(length(records))
  gradient = matrix(0, length(records), 2 * p)

  # Keep each configuration's complete coefficient vector in the calculation.
  for (i in seq_along(records)) {
    record = records[[i]]
    row = if (coefficient) diag(p)[record$parameter, ] else problem$X[record$configuration, ]
    mean = sum(row * parameters$mean)
    variance = sum((row * parameters$sd)^2)
    sd = sqrt(variance)
    mean_gradient = row * problem$sd
    if (all(row == 0)) {
      value[i] = prior_joint_constant_value(record)
      next
    }
    if (record$type %in% c("mean", "mode")) {
      value[i] = mean
      gradient[i, seq_len(p)] = mean_gradient
    } else if (record$type == "sd") {
      value[i] = sd
      gradient[i, p + seq_len(p)] = (row * parameters$sd)^2 / sd
    } else {
      z = (record$threshold - mean) / sd
      sign = if (record$operator %in% c("<", "<=")) 1 else -1
      value[i] = if (sign == 1) stats::pnorm(z) else stats::pnorm(z, lower.tail = FALSE)
      gradient[i, seq_len(p)] = -sign * stats::dnorm(z) * mean_gradient / sd
      gradient[i, p + seq_len(p)] = -sign * stats::dnorm(z) * z * (row * parameters$sd)^2 / variance
    }
  }
  list(value = value, gradient = gradient)
}

prior_joint_kl = function(theta, problem) {
  p = length(problem$mean)
  u = theta[seq_len(p)]
  l = theta[p + seq_len(p)]
  list(value = sum(u^2 + expm1(2 * l) - 2 * l) / 2, gradient = c(u, expm1(2 * l)))
}

# Belief sacrifices are distances outside acceptance bands in dimensionless
# coordinates. A probability's probit is its threshold in Normal SD units;
# use it directly instead of inverting a rounded/saturated CDF value.
prior_joint_coefficient_metric = function(theta, problem) {
  records = problem$coefficient
  raw = prior_joint_values(theta, problem, TRUE)
  parameters = prior_joint_parameters(theta, problem)
  p = length(problem$mean)
  value = lower = upper = numeric(length(records))
  gradient = raw$gradient
  for (i in seq_along(records)) {
    record = records[[i]]
    j = record$parameter
    if (record$type != "probability") {
      value[i] = (raw$value[i] - record$target) / problem$sd[j]
      lower[i] = -record$tolerance / problem$sd[j]
      upper[i] = -lower[i]
      gradient[i, ] = gradient[i, ] / problem$sd[j]
      next
    }
    sign = if (record$operator %in% c("<", "<=")) 1 else -1
    z = sign * (record$threshold - parameters$mean[j]) / parameters$sd[j]
    target = stats::qnorm(record$target)
    value[i] = z - target
    lower[i] = stats::qnorm(max(0, record$target - record$tolerance)) - target
    upper[i] = stats::qnorm(min(1, record$target + record$tolerance)) - target
    gradient[i, ] = 0
    gradient[i, j] = -sign * problem$sd[j] / parameters$sd[j]
    gradient[i, p + j] = -z
  }
  list(value = value, lower = lower, upper = upper, gradient = gradient,
    violation = pmax(0, lower - value, value - upper))
}

# Scalar objective and inequalities (g <= 0), with analytic derivatives.
prior_joint_system = function(theta, problem, tier, coefficient_limit = 0, coefficient_scale = 1) {
  p = length(problem$mean) * 2L
  parameters = theta[seq_len(p)]
  outcome = prior_joint_values(parameters, problem)
  target = vapply(problem$outcome, `[[`, numeric(1), "target")
  tolerance = vapply(problem$outcome, `[[`, numeric(1), "tolerance")
  difference = outcome$value - target

  # A small interior margin protects final checks from floating-point rounding.
  # Retry without it when the tolerance boxes have only a narrow intersection.
  width = tolerance * (1 - problem$interior)
  constant = vapply(problem$outcome, function(r) all(problem$X[r$configuration, ] == 0), logical(1))
  width[constant] = tolerance[constant]
  g = c((difference - width) / problem$units, (-difference - width) / problem$units)
  jacobian = rbind(outcome$gradient / problem$units, -outcome$gradient / problem$units)
  if (length(problem$coefficient) && tier != 1L) {
    coefficient = prior_joint_coefficient_metric(parameters, problem)
    upper = which(is.finite(coefficient$upper))
    lower = which(is.finite(coefficient$lower))
    allowance = if (tier == 2L) coefficient_scale * theta[p + 1L] else coefficient_limit
    g = c(g, coefficient$value[upper] - coefficient$upper[upper] - allowance,
      coefficient$lower[lower] - coefficient$value[lower] - allowance)
    jacobian = rbind(jacobian, coefficient$gradient[upper, , drop = FALSE],
      -coefficient$gradient[lower, , drop = FALSE])
  }
  if (tier == 2L) {
    jacobian = cbind(jacobian, c(rep(0, 2 * length(target)),
      rep(-coefficient_scale, length(upper) + length(lower))))
    g = c(g, -theta[p + 1L])
    jacobian = rbind(jacobian, c(rep(0, p), -1))
    return(list(value = theta[p + 1L], gradient = c(rep(0, p), 1), g = g, jacobian = jacobian))
  }
  kl = prior_joint_kl(parameters, problem)
  list(value = kl$value, gradient = kl$gradient, g = g, jacobian = jacobian)
}

prior_joint_augmented = function(theta, problem, tier, coefficient_limit, coefficient_scale, lambda, rho) {
  system = prior_joint_system(theta, problem, tier, coefficient_limit, coefficient_scale)
  if (any(!is.finite(c(system$value, system$gradient, system$g, system$jacobian)))) return(list(value = 1e100, gradient = rep(0, length(theta))))
  multiplier = pmax(0, lambda + rho * system$g)
  list(value = system$value + sum(multiplier^2 - lambda^2) / (2 * rho),
    gradient = system$gradient + as.vector(crossprod(system$jacobian, multiplier)))
}
prior_joint_augmented_value = function(theta, ...) prior_joint_augmented(theta, ...)$value
prior_joint_augmented_gradient = function(theta, ...) prior_joint_augmented(theta, ...)$gradient

prior_joint_solve = function(start, problem, tier = 3L, coefficient_limit = 0, coefficient_scale = 1) {
  reference_start = start
  lambda = rep(0, length(prior_joint_system(start, problem, tier, coefficient_limit, coefficient_scale)$g))
  rho = 10
  previous = Inf
  for (iteration in seq_len(24L)) {
    fit = stats::optim(
      start,
      prior_joint_augmented_value,
      prior_joint_augmented_gradient,
      problem = problem,
      tier = tier,
      coefficient_limit = coefficient_limit,
      coefficient_scale = coefficient_scale,
      lambda = lambda,
      rho = rho,
      method = "BFGS",
      control = list(maxit = 600, reltol = 1e-13)
    )
    current = prior_joint_system(fit$par, problem, tier, coefficient_limit, coefficient_scale)
    violation = max(c(0, current$g))
    step = max(abs(fit$par - start))
    lambda = pmax(0, lambda + rho * current$g)
    start = fit$par
    if (violation < 1e-12 && step < 1e-7 && fit$convergence == 0L) break
    if (violation > previous * .3) rho = min(rho * 10, 1e14)
    previous = violation
  }

  # Restore numerical feasibility after the augmented-Lagrangian iterations.
  # Baseline-scaled constraints can be tiny in dollar-valued examples; an
  # optimizer's absolute stopping rule is insufficient for belief tolerances.
  for (polish in seq_len(20L)) {
    current = prior_joint_system(start, problem, tier, coefficient_limit, coefficient_scale)

    # Outcome acceptance can be far tighter than generic optimizer accuracy.
    # Do not let roundoff-size coefficient-bound errors repeatedly undo an
    # outcome correction. This is numerical allowance, not a belief weight;
    # zero coefficient sacrifice is still checked against raw tolerances.
    allowance = rep(0, length(current$g))
    if ((tier == 2L && tail(start, 1) > 1e-8) || (tier == 3L && coefficient_limit > 0)) {
      coefficient_rows = seq_along(current$g) > 2 * length(problem$outcome)
      allowance[coefficient_rows] = 1e-10
    }
    active = which(current$g > allowance)
    if (!length(active)) break
    J = current$jacobian[active, , drop = FALSE]
    singular = svd(J)
    keep = singular$d > max(singular$d) * 1e-12
    if (!any(keep)) break
    correction = as.vector(singular$v[, keep, drop = FALSE] %*%
      (crossprod(singular$u[, keep, drop = FALSE], current$g[active]) / singular$d[keep]))
    start = start - correction
  }
  current = prior_joint_system(start, problem, tier, coefficient_limit, coefficient_scale)
  result = list(
    par = start,
    convergence = fit$convergence,
    iterations = iteration,
    violation = max(c(0, current$g)),
    objective = current$value,
    stationarity = prior_joint_stationarity(current)
  )
  if (tier != 3L || prior_joint_stage_verified(result)) return(result)

  # The penalty solve can stall on a very narrow near-minimax face. Refine
  # equivalent constraints while preserving the priorities and acceptance checks.
  refined = tryCatch(prior_joint_face_refine(reference_start, problem, coefficient_limit, coefficient_scale), error = function(error) NULL)
  if (is.null(refined)) result else refined
}

prior_joint_stage_verified = function(stage) {
  stage$convergence == 0 && is.finite(stage$objective) && stage$violation < 1e-8 &&
    is.finite(stage$stationarity) && stage$stationarity < 1e-4
}

# A nonnegative multiplier certificate checks local first-order optimality.
# Dropping negative multipliers is conservative: failure leaves the problem
# unresolved rather than accepting a merely feasible, nonstationary solution.
prior_joint_stationarity = function(system) {
  active = which(system$g >= -1e-7)
  residual = system$gradient
  while (length(active)) {
    J = system$jacobian[active, , drop = FALSE]
    singular = svd(t(J))
    keep = singular$d > max(singular$d) * 1e-12
    if (!any(keep)) break
    multiplier = as.vector(singular$v[, keep, drop = FALSE] %*%
      (crossprod(singular$u[, keep, drop = FALSE], -system$gradient) / singular$d[keep]))
    if (all(multiplier >= 0)) {
      residual = c(system$gradient + as.vector(crossprod(J, multiplier)),
        multiplier * system$g[active])
      break
    }
    active = active[-which.min(multiplier)]
  }
  max(abs(residual))
}

prior_joint_projection = function(X, target, mean, sd) {
  A = sweep(X, 2, sd, "*")
  singular = svd(A)
  keep = singular$d > max(singular$d) * 1e-12
  adjustment = if (any(keep)) as.vector(singular$v[, keep, drop = FALSE] %*%
    (crossprod(singular$u[, keep, drop = FALSE], target - as.vector(X %*% mean)) / singular$d[keep])) else rep(0, length(mean))
  mean + sd * adjustment
}

prior_joint_problem = function(object) {
  X = prior_joint_design(object, object$outcomes)
  names = colnames(X)
  baseline = object$selected
  if (any(vapply(baseline[names], is.null, logical(1)))) stop("Resolve coefficient-only priors before adding outcome beliefs.")
  if (any(vapply(baseline[names], function(p) p$distribution != "normal", logical(1)))) {
    stop("Joint outcome elicitation requires independent Normal priors, including the centered intercept; Student's t priors are unsupported.")
  }
  mean = vapply(baseline[names], `[[`, numeric(1), "mean")
  sd = vapply(baseline[names], `[[`, numeric(1), "sd")
  outcome = list()
  coefficient = list()
  for (i in seq_along(object$outcomes)) for (record in object$outcomes[[i]]$beliefs) {
    record$configuration = i
    outcome = c(outcome, list(record))
  }
  for (i in seq_along(names)) for (record in object$beliefs[[names[i]]]) {
    record$parameter = i
    coefficient = c(coefficient, list(record))
  }
  variance = as.vector(X^2 %*% sd^2)
  units = vapply(outcome, function(r) if (r$type == "probability") 1 else if (variance[r$configuration] == 0) 1 else sqrt(variance[r$configuration]), numeric(1))
  coefficient_units = vapply(coefficient, function(r) if (r$type == "probability") 1 else sd[r$parameter], numeric(1))
  problem = list(
    X = X,
    mean = mean,
    sd = sd,
    outcome = outcome,
    coefficient = coefficient,
    units = units,
    coefficient_units = coefficient_units,
    baseline = baseline,
    interior = .001
  )
  if (any(vapply(object$outcomes, function(specification) !is.null(specification$difference), logical(1)))) problem$has_difference = TRUE
  problem
}

prior_joint_basic_compatibility = function(problem) {
  for (i in seq_len(nrow(problem$X))) {
    if (all(problem$X[i, ] == 0)) {
      indices = which(vapply(problem$outcome, function(r) r$configuration == i, logical(1)))
      bad = indices[vapply(problem$outcome[indices], function(r)
        abs(prior_joint_constant_value(r) - r$target) > r$tolerance, logical(1))]
      if (length(bad)) return(list(issue = "zero contrast",
        message = "Difference beliefs are incompatible: these configurations have an expected difference of exactly zero.",
        configurations = i, records = bad,
        explanation = "The model-matrix contrast is zero. The expected difference has no uncertainty; residual noise is not included.",
        technical = list(method = "Deterministic zero contrast", contrast = problem$X[i, ])))
      next
    }
    rows = which(apply(problem$X, 1, function(x) identical(unname(x), unname(problem$X[i, ]))))
    records = Filter(function(r) r$configuration %in% rows, problem$outcome)
    check = prior_check_probability_constraints(records)
    if (check$contradictory) {
      return(list(
        issue = "probability",
        message = check$message,
        configurations = rows,
        records = which(vapply(problem$outcome, function(r)
          r$configuration %in% rows && r$type == "probability", logical(1))),
        explanation = "The stated probabilities conflict, regardless of the selected coefficient distributions.",
        technical = list(method = "CDF probability consistency")
      ))
    }
    location = Filter(function(r) r$type %in% c("mean", "mode"), records)
    if (length(location) > 1L) {
      lower = max(vapply(location, function(r) r$target - r$tolerance, numeric(1)))
      upper = min(vapply(location, function(r) r$target + r$tolerance, numeric(1)))
      unit = sqrt(sum((problem$X[i, ] * problem$sd)^2))
      if (lower - upper > 1e-12 * unit) {
        return(list(
          issue = "location",
          message = "Outcome location constraints are incompatible at the same predictor configuration.",
          configurations = rows,
          records = which(vapply(problem$outcome, function(r)
            r$configuration %in% rows && r$type %in% c("mean", "mode"), logical(1))),
          explanation = "The requested expected-value intervals at the same predictor values do not overlap.",
          technical = list(method = "Disjoint location acceptance intervals", lower = lower, upper = upper)
        ))
      }
    }
  }
  location = Filter(function(r) r$type %in% c("mean", "mode") &&
    any(problem$X[r$configuration, ] != 0), problem$outcome)
  if (!length(location)) return(NULL)
  A = problem$X[vapply(location, `[[`, integer(1), "configuration"), , drop = FALSE]
  b = vapply(location, `[[`, numeric(1), "target")
  tol = vapply(location, `[[`, numeric(1), "tolerance")

  # Test the certificate in dimensionless coordinates as well: neither raw
  # predictor units nor the outcome currency should change a conflict label.
  scaled = sweep(A, 2, problem$sd, "*")
  units = sqrt(rowSums(scaled^2))
  b = (b - as.vector(A %*% problem$mean)) / units
  A = scaled / units
  tol = tol / units
  A = A[, colSums(abs(A)) > 0, drop = FALSE]

  # Do not infer impossibility by truncating a tiny but nonzero singular
  # direction: an ill-conditioned mean problem may still have a large solution.
  if (nrow(A) <= ncol(A)) return(NULL)
  singular = svd(A, nu = 0, nv = 0)$d
  if (min(singular) <= sqrt(.Machine$double.eps) * max(singular)) return(NULL)
  projected = prior_joint_projection(A, b, rep(0, ncol(A)), rep(1, ncol(A)))
  residual = b - as.vector(A %*% projected)

  # A left-null-space certificate proves that the tolerance boxes cannot meet.
  if (max(abs(crossprod(A, residual))) < 1e-8 * max(1, sum(abs(A) * max(abs(residual)))) &&
      abs(sum(residual * b)) > sum(abs(residual) * tol) + 1e-9 * max(1, sum(abs(residual * b)))) {
    return(list(
      issue = "location",
      message = "Outcome location constraints are incompatible with the model matrix, even within their tolerances.",
      configurations = unique(vapply(location, `[[`, integer(1), "configuration")),
      records = which(vapply(problem$outcome, function(r) r$type %in% c("mean", "mode"), logical(1))),
      explanation = "The requested expected values cannot all hold in this regression, even within their acceptance intervals.",
      technical = list(method = "Location tolerance-box separation", design = A, residual = residual)
    ))
  }
  NULL
}

prior_joint_details = function(problem, theta, coefficient = FALSE) {
  records = if (coefficient) problem$coefficient else problem$outcome
  values = prior_joint_values(theta, problem, coefficient)$value
  target = vapply(records, `[[`, numeric(1), "target")
  tolerance = vapply(records, `[[`, numeric(1), "tolerance")
  residual = (values - target) / tolerance
  discrepancy = values - target

  # Converted standardized constraints retain the researcher's entered units.
  factors = vapply(records, function(r) if (coefficient && r$type != "probability" &&
    !is.null(r$standardization_factor)) r$standardization_factor else 1, numeric(1))
  details = data.frame(
    id = vapply(records, `[[`, numeric(1), "id"),
    target = target / factors,
    achieved = values / factors,
    discrepancy = discrepancy / factors,
    tolerance = tolerance / factors,
    residual = residual,
    satisfied = is.finite(residual) & abs(residual) <= 1
  )
  if (!coefficient) return(details)
  metric = prior_joint_coefficient_metric(theta, problem)
  details$metric = vapply(records, function(r) if (r$type == "probability") "probit" else "baseline SD", character(1))
  details$normalized_discrepancy = metric$value
  details$normalized_violation = metric$violation
  details
}

# Boundary probability targets use an interior starting point of the declared
# acceptance band. Verification still uses the original target and tolerance.
prior_joint_start = function(problem) {
  problem = prior_joint_ordered_problem(problem)
  targets = vapply(problem$outcome, function(r) {
    if (all(problem$X[r$configuration, ] == 0)) return(0)
    if (r$type != "probability") return(r$target)
    probability = r$target
    if (probability <= 0 || probability >= 1) {
      probability =
        (max(0, r$target - r$tolerance) + min(1, r$target + r$tolerance)) / 2
    }
    z = stats::qnorm(probability, lower.tail = r$operator %in% c("<", "<="))
    if (!is.finite(z)) stop("Probability acceptance band cannot be initialized at finite numerical precision.")
    r$threshold - z * sqrt(sum((problem$X[r$configuration, ] * problem$sd)^2))
  }, numeric(1))
  A = problem$X[vapply(problem$outcome, `[[`, integer(1), "configuration"), , drop = FALSE]
  location = prior_joint_projection(A, targets, problem$mean, problem$sd)
  c((location - problem$mean) / problem$sd, rep(0, length(problem$sd)))
}

prior_joint_attempt = function(start, problem, ...) {
  problem = prior_joint_ordered_problem(problem)
  tryCatch(prior_joint_solve(start, problem, ...), error = function(error)
    list(par = start, convergence = 1L, iterations = 0L, violation = Inf, objective = Inf,
      stationarity = Inf, message = conditionMessage(error)))
}

prior_joint_reconcile = function(object) {
  object$joint = NULL
  problem = prior_joint_problem(object)

  # Only a verified mathematical certificate establishes incompatibility.
  assessment = tryCatch(prior_joint_compatibility(problem), error = identity)
  certificate = if (inherits(assessment, "error")) NULL else assessment
  conflict = if (is.null(certificate)) NULL else certificate$message
  base = rep(0, 2 * length(problem$mean))
  joint = list(
    baseline = problem$baseline,
    beliefs = object$beliefs,
    specifications = object$outcomes,
    design = problem$X,
    problem = problem,
    status = "unresolved",
    message = conflict,
    compatibility = certificate,
    compatibility_note = if (inherits(assessment, "error")) conditionMessage(assessment) else NULL
  )
  if (!is.null(conflict)) {
    joint$status = "incompatible"
    object$joint = joint
    warning(conflict, call. = FALSE)
    return(object)
  }

  # A projection provides a deterministic, unit-aware starting point, including
  # single probability constraints that do not identify coefficient parameters.
  start = tryCatch(prior_joint_start(problem), error = identity)
  if (inherits(start, "error")) {
    joint$message = paste("Numerical initialization failed; no mathematical incompatibility was established.", conditionMessage(start))
    object$joint = joint
    warning(joint$message, call. = FALSE)
    return(object)
  }

  initial = prior_joint_attempt(start, problem, tier = 1L)
  if (!all(prior_joint_details(problem, initial$par)$satisfied) || initial$violation > 1e-8) {
    problem$interior = 0
    initial = prior_joint_attempt(start, problem, tier = 1L)
  }
  joint$problem = problem
  stages = list(outcome = initial)
  theta = initial$par
  outcome = prior_joint_details(problem, theta)
  if (all(outcome$satisfied) && length(problem$coefficient)) {
    scale = max(1, prior_joint_coefficient_metric(theta, problem)$violation)
    second = prior_joint_attempt(c(theta, 1), problem, tier = 2L, coefficient_scale = scale)
    theta2 = head(second$par, length(base))
    limit = max(c(0, prior_joint_coefficient_metric(theta2, problem)$violation))

    # A relative numerical allowance opens the nearly optimal minimax face
    # for KL optimization. Zero sacrifice stays exactly zero.
    third = prior_joint_attempt(theta2, problem, tier = 3L, coefficient_limit = limit * (1 + 1e-6), coefficient_scale = scale)
    if (!is.null(third$minimax)) {
      second = third$minimax
      theta2 = head(second$par, length(base))
      limit = third$coefficient_limit
    }
    joint$coefficient_limit = limit
    stages = c(stages, list(coefficient = second, KL = third))
    theta = third$par
  }
  joint$outcome = prior_joint_details(problem, theta)
  joint$coefficient = prior_joint_details(problem, theta, TRUE)
  joint$stages = stages
  joint$theta = theta
  joint$kl = prior_joint_kl(theta, problem)$value
  joint$reference_used = TRUE
  joint$underdetermined = qr(prior_joint_values(theta, problem)$gradient / problem$units)$rank < length(theta)
  verified = all(joint$outcome$satisfied) && all(vapply(stages, prior_joint_stage_verified, logical(1)))
  if (!is.null(joint$coefficient_limit)) {
    achieved = max(c(0, joint$coefficient$normalized_violation))
    bound = joint$coefficient_limit * (1 + 1e-6)
    tier_order = achieved <= bound + 1e-10 &&
      (bound > 0 || all(joint$coefficient$satisfied)) &&
      joint$kl <= prior_joint_kl(theta2, problem)$value + 1e-8 * max(1, prior_joint_kl(theta2, problem)$value)
    verified = verified && tier_order
  }

  # Recheck the actual distributions before treating the result as selected.
  parameters = prior_joint_parameters(theta, problem)
  verified = verified && all(is.finite(c(parameters$mean, parameters$sd))) && all(parameters$sd > 0)
  if (!verified) {
    joint$message = "The numerical solver could not verify a feasible joint prior. No mathematical incompatibility was established; inspect diagnostics or revisit the model and acceptance tolerances."
    object$joint = joint
    warning(joint$message, call. = FALSE)
    return(object)
  }
  for (i in seq_along(problem$mean)) object$selected[[names(problem$mean)[i]]] =
    list(distribution = "normal", mean = unname(parameters$mean[i]), sd = unname(parameters$sd[i]))
  joint$adjusted = object$selected
  joint$status = if (all(joint$coefficient$satisfied)) "compatible" else "reconciled"
  object$joint = joint
  if (joint$status == "reconciled") warning("Some coefficient beliefs were relaxed to satisfy outcome beliefs; inspect summary(..., diagnostics=TRUE).", call. = FALSE)
  object
}

# Recompute achieved quantities; do not trust stale reconciliation annotations.
prior_joint_validate = function(object) {
  joint = object$joint
  if (is.null(joint)) return(invisible(object))
  if (!joint$status %in% c("compatible", "reconciled")) stop("Outcome reconciliation is ", joint$status, ": ", joint$message)
  if (!identical(object$beliefs, joint$beliefs) || !identical(object$outcomes, joint$specifications) ||
      !isTRUE(all.equal(object$selected, joint$adjusted))) {
    stop("Joint prior metadata or selected distributions changed; re-elicit the priors.")
  }
  prior_validate_standardization_data(object)
  baseline = object
  baseline$selected = joint$baseline
  rebuilt = prior_joint_problem(baseline)
  fields = setdiff(names(rebuilt), "interior")
  if (!isTRUE(all.equal(rebuilt[fields], joint$problem[fields]))) stop("Reconciliation baseline or constraints changed; re-elicit the priors.")
  parameters = prior_joint_parameters(joint$theta, joint$problem)
  for (i in seq_along(parameters$mean)) {
    selected = object$selected[[names(parameters$mean)[i]]]
    if (!isTRUE(all.equal(unname(c(selected$mean, selected$sd)),
        unname(c(parameters$mean[i], parameters$sd[i]))))) {
      stop("Adjusted priors do not match the reconciliation parameters.")
    }
  }
  X = prior_prediction_matrix(object, object$model$data)
  if (!isTRUE(all.equal(X, object$model$X))) stop("Model data changed since outcome elicitation; identify priors again.")
  if (!isTRUE(all.equal(prior_joint_design(object, object$outcomes), joint$design))) stop("Outcome prediction design changed; re-elicit the priors.")
  outcome = prior_joint_details(joint$problem, joint$theta)
  coefficient = prior_joint_details(joint$problem, joint$theta, TRUE)
  if (!all(outcome$satisfied) || !isTRUE(all.equal(outcome, joint$outcome)) ||
      !isTRUE(all.equal(coefficient, joint$coefficient))) {
    stop("Reconciliation constraints or diagnostics do not validate.")
  }
  invisible(object)
}

prior_joint_attempted_priors = function(joint) {
  if (is.null(joint$problem) || length(joint$theta) != 2 * length(joint$problem$mean)) return(list())
  parameters = prior_joint_parameters(joint$theta, joint$problem)
  if (any(!is.finite(c(parameters$mean, parameters$sd))) || any(parameters$sd <= 0)) return(list())
  result = lapply(seq_along(parameters$mean), function(i) list(distribution = "normal",
    mean = unname(parameters$mean[i]), sd = unname(parameters$sd[i])))
  stats::setNames(result, names(joint$problem$mean))
}

prior_joint_diagnostics = function(x, verified = TRUE) {
  object = x$object
  joint = object$joint
  cat("\nTechnical diagnostics (", joint$status, ")\n", sep = "")
  if (!verified) prior_console_line("Unverified reconciliation: printed available priors are retained fallbacks; attempted evaluations below are not selected priors.")
  print(object)
  if (!is.null(joint$message)) prior_console_line(joint$message)
  if (!is.null(joint$compatibility_note)) prior_console_line(paste("Compatibility check unavailable:", joint$compatibility_note))
  prior_joint_summary_conflict(object, TRUE)
  if (!verified) cat("\nAttempted outcome evaluations (not selected):\n")
  index = 0L
  for (configuration in seq_along(object$outcomes)) {
    specification = object$outcomes[[configuration]]
    prior_joint_summary_configuration(specification)
    if (!is.null(specification$difference)) {
      cat("Model-matrix contrast (to - from):\n")
      print(joint$design[configuration, , drop = FALSE])
    }
    for (record in specification$beliefs) {
      index = index + 1L
      prior_console_line(record$label)
      if (!is.null(joint$outcome)) print(joint$outcome[index, , drop = FALSE], row.names = FALSE, digits = 8)
    }
  }
  cat("\nCoefficient reconciliation (centered-intercept parameterization):\n")
  cat("Discrepancy/tolerance: entered units; residual: discrepancy / numerical tolerance.\n")
  cat("Comparison metric: baseline-SD distances for locations/SDs, probit distances for probabilities.\n")
  cat("Normalized violation: distance outside the tolerance band; zero does not override the raw satisfaction check.\n")
  attempted = prior_joint_attempted_priors(joint)
  if (verified) cat("Before -> after: original coefficient prior -> verified reconciled prior.\n")
  if (!verified) cat("Attempted coefficient belief evaluations are unverified and refer to the attempted priors, not the retained fallbacks.\n")
  for (name in names(joint$problem$mean)) {
    before = joint$baseline[[name]]
    after = if (verified) object$selected[[name]] else attempted[[name]]
    if (verified) cat(name, ":", prior_format(before), "->", prior_format(after), "\n") else {
      cat(name, ":\n")
      cat("  Original coefficient prior:", prior_format(before), "\n")
      cat("  Retained fallback prior:", prior_format(object$selected[[name]]), "\n")
      if (!is.null(after)) cat("  Attempted numerical prior (not selected):", prior_format(after), "\n")
    }
    if (!is.null(after)) cat(if (verified) "  Change in mean:" else "  Attempted change in mean:",
      after$mean - before$mean, "; change in SD:", after$sd - before$sd, "\n")
    records = object$beliefs[[name]]
    for (record in records) {
      row = if (is.null(joint$coefficient)) data.frame() else
        joint$coefficient[joint$coefficient$id == record$id, , drop = FALSE]
      prior_console_line(paste(record$label, if (!verified) "[attempted; not selected]" else
        if (nrow(row) && row$satisfied) "[satisfied]" else "[relaxed / unresolved]"))
      if (nrow(row)) print(row, row.names = FALSE, digits = 8)
    }
  }
  if (!is.null(joint$coefficient_limit)) cat("\nTier-2 maximum dimensionless violation:", joint$coefficient_limit, "\n")
  if (!is.null(joint$kl)) cat(if (verified) "\nKL(adjusted || baseline):" else "\nKL(attempted || baseline):", joint$kl, "\n")
  cat("Reference anchoring: coefficient-only selections, with backend defaults where not elicited.\n")
  cat("Underdetermined constraints use the stated minimax/KL criterion; numerical search does not certify a unique global optimum.\n")
  for (stage in names(joint$stages)) if (!is.null(joint$stages[[stage]]$message))
    prior_console_line(paste(stage, ":", joint$stages[[stage]]$message))
  for (stage in names(joint$stages)) if (!is.null(joint$stages[[stage]]$method))
    prior_console_line(paste(stage, "method:", joint$stages[[stage]]$method))
  if (length(joint$stages)) {
    print(do.call(rbind, lapply(joint$stages, function(s) data.frame(
      convergence = s$convergence, iterations = s$iterations, scaled_violation = s$violation,
      stationarity = s$stationarity, objective = s$objective))))
  }
  invisible(x)
}

prior_joint_summary_number = function(value, probability = FALSE, mark = ",") {
  if (length(value) != 1L || !is.finite(value)) return("Not available")
  if (probability) value = 100 * value
  scientific = abs(value) > 1e12 || (value != 0 && abs(value) < 1e-4)
  result = format(signif(value, if (probability) 4 else 6), trim = TRUE, big.mark = mark, scientific = scientific)
  paste0(result, if (probability) "%" else "")
}

prior_joint_summary_belief = function(record, outcome = "value") {
  if (record$type == "probability") {
    threshold = if (is.null(record$entered)) record$threshold else record$entered$threshold
    label = paste0("P(", outcome, " ", record$operator, " ", prior_joint_summary_number(threshold), ")")
  } else label = switch(record$type, mode = paste("Most likely", outcome),
    mean = paste("Mean", outcome), sd = paste("SD of", outcome), record$label)
  if (identical(record$scale, "standardized")) label = paste0(label, " [standardized]")
  label
}

prior_joint_summary_distribution = function(prior) {
  if (is.null(prior)) return("Not selected")
  if (prior$distribution != "normal") return(prior_format(prior))
  paste0("N(", prior_joint_summary_number(prior$mean, mark = ""), ",",
    prior_joint_summary_number(prior$sd, mark = ""), ")")
}

prior_joint_summary_outcomes = function(object, ready) {
  outcome = paste(deparse(object$model$formula[[2]]), collapse = " ")
  cat("\nOutcome:", outcome, "\n")
  index = 0L
  for (specification in object$outcomes) {
    if (!is.null(specification$difference)) cat("\nDifference in expected", outcome, "(to - from)\n")
    prior_joint_summary_configuration(specification)
    rows = list()
    for (record in specification$beliefs) {
      index = index + 1L
      achieved = if (ready) prior_joint_summary_number(object$joint$outcome$achieved[index],
        probability = record$type == "probability") else "Not verified"
      label = prior_joint_summary_belief(record, if (!is.null(specification$difference)) "expected difference" else paste("expected", outcome))
      if (!is.null(object$joint$compatibility)) label = paste0(length(rows) + 1L, ". ", label)
      rows[[length(rows) + 1L]] = data.frame(
        Belief = label,
        Requested = prior_joint_summary_number(record$target, probability = record$type == "probability"),
        Achieved = achieved)
    }
    cat("\nOutcome beliefs\n")
    print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
    cat("\n")
  }
  satisfied = if (ready) "Yes" else if (object$joint$status == "incompatible") "No (conflicting beliefs)" else "Not verified"
  cat("Outcome beliefs satisfied:", satisfied, "\n")
}

# Build one display group per changed coefficient, followed only by beliefs
# whose recorded satisfaction checks failed. No belief evaluation is repeated.
prior_joint_summary_comparison = function(object) {
  joint = object$joint
  parameters = names(joint$problem$mean)
  changed = parameters[!vapply(parameters, function(name)
    identical(object$selected[[name]], joint$baseline[[name]]), logical(1))]
  rows = list()
  for (parameter in changed) {
    if (length(rows)) rows[[length(rows) + 1L]] = c("", "", "")
    rows[[length(rows) + 1L]] = c(parameter, "", "")
    rows[[length(rows) + 1L]] = c("  Prior", prior_joint_summary_distribution(joint$baseline[[parameter]]),
      prior_joint_summary_distribution(object$selected[[parameter]]))
    for (record in object$beliefs[[parameter]]) {
      achieved = joint$coefficient[joint$coefficient$id == record$id, , drop = FALSE]
      if (!nrow(achieved) || achieved$satisfied) next
      rows[[length(rows) + 1L]] = c(paste0("  ", prior_joint_summary_belief(record)),
        prior_joint_summary_number(achieved$target, probability = record$type == "probability"),
        prior_joint_summary_number(achieved$achieved, probability = record$type == "probability"))
    }
  }
  if (!length(rows)) return(NULL)
  result = as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(result) = c("Parameter / Belief", "Before", "After")
  result
}

# Wrap cells, including unusually long parameter names, without truncating data.
prior_joint_summary_cell = function(text, width) {
  if (!nzchar(text) || nchar(text) <= width) return(text)
  prefix = if (startsWith(text, "  ")) "  " else ""
  width = max(1, width - nchar(prefix))
  words = strwrap(trimws(text), width = width)
  if (!length(words)) return("")
  parts = lapply(words, function(word) {
    starts = seq.int(1, nchar(word), by = width)
    substring(word, starts, pmin(nchar(word), starts + width - 1))
  })
  paste0(prefix, unlist(parts, use.names = FALSE))
}

prior_joint_summary_row = function(row, widths) {
  if (all(!nzchar(unlist(row)))) { cat("\n")
    return(invisible(NULL)) }
  cells = Map(prior_joint_summary_cell, row, widths)
  for (i in seq_len(max(lengths(cells)))) {
    line = vapply(cells, function(cell) if (length(cell) >= i) cell[i] else "", character(1))
    line = paste(sprintf(paste0("%-", widths, "s"), line), collapse = "  ")
    cat(sub(" +$", "", line), "\n", sep = "")
  }
}

prior_joint_summary_coefficients = function(object, ready) {
  cat("\nCoefficient priors and beliefs\n")
  if (!ready) {
    cat("No verified reconciliation to report.\n")
    return(invisible(NULL))
  }
  rows = prior_joint_summary_comparison(object)
  if (is.null(rows)) { cat("No coefficient priors changed.\n")
    return(invisible(NULL)) }
  widths = pmax(nchar(names(rows)), vapply(rows, function(column) max(nchar(column)), integer(1)))
  console = max(40L, getOption("width", 80L))
  widths[1] = min(widths[1], max(12, console - sum(widths[-1]) - 4))
  if (sum(widths) + 4 > console) {
    widths[2] = min(widths[2], floor((console - widths[1] - 4) / 2))
    widths[3] = min(widths[3], console - widths[1] - widths[2] - 4)
  }
  prior_joint_summary_row(as.list(names(rows)), widths)
  cat(strrep("-", sum(widths) + 4), "\n", sep = "")
  for (i in seq_len(nrow(rows))) prior_joint_summary_row(as.list(rows[i, ]), widths)
  cat("Other priors unchanged.\n")
}

prior_joint_summary = function(x) {
  object = x$object
  error = tryCatch({ prior_predictive_validate(object)
    NULL }, error = conditionMessage)
  ready = is.null(error)
  status = if (ready) "Successful" else switch(object$joint$status,
    incompatible = "Incompatible", unresolved = "Not resolved", "Needs review")
  cat("Prior reconciliation:", status, "\n")
  if (object$backend$package == "brms") prior_console_line("brms backend; starting priors are proper flexbayes defaults derived from rstanarm, overriding native brms defaults.")
  prior_joint_summary_outcomes(object, ready)
  prior_joint_summary_coefficients(object, ready)
  prior_r_squared_summary(object, x$diagnostics)
  if (isTRUE(object$joint$retained)) {
    cat("\nAvailable priors (retained; outcome beliefs not reconciled):\n")
    print(object)
  }
  if (length(object$beliefs$sigma)) {
    cat("\nResidual-scale prior:", prior_format(object$selected$sigma), "\n")
    for (i in seq_along(object$beliefs$sigma)) {
      prior_console_line(paste0("  ", i, ". ", object$beliefs$sigma[[i]]$label))
    }
    if (!is.null(object$candidates$sigma$message)) prior_console_line(object$candidates$sigma$message)
    prior_scale_summary(object, x$diagnostics)
  }
  if (!ready) {
    if (object$joint$status == "incompatible") prior_joint_summary_conflict(object)
    reason = if (identical(object$joint$attempted, FALSE)) object$joint$message else if (object$joint$status == "unresolved")
      "The solver could not verify a feasible prior; this does not prove your answers are impossible." else if (object$joint$status == "incompatible")
      "Review the conflicting answers or the model's restrictions before using these priors." else error
    prior_console_line(paste("Not ready to use:", reason))
  }
  if (x$diagnostics) prior_joint_diagnostics(x, ready) else
    prior_console_line("For mathematical details: summary(priors, diagnostics = TRUE).")
  invisible(x)
}
