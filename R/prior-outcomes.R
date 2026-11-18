prior_outcome_scale = function(scale) !is.null(scale) && scale %in% c("probability", "response", "odds_ratio", "rate_ratio")

prior_outcome_validate = function(record, object, parameter) {
  if (prior_ratio_scale(record$scale)) return(prior_ratio_validate(record, object, parameter))
  if (parameter != "(Intercept)") stop("Outcome-scale beliefs currently support the centered intercept only.")
  expected = switch(object$model$family$family, binomial = "probability", poisson = "response", "unsupported")
  if (record$scale != expected) stop("This family does not support the requested outcome scale: ", record$scale)

  # Targets must lie within the outcome domain before transformation.
  value = if (record$type == "probability") record$threshold else record$target
  if (record$scale == "probability") {
    if (record$type == "sd") {
      if (value >= .5) stop("A probability SD must be strictly below 0.5.")
    } else if (record$type == "probability") {
      if (value < 0 || value > 1) stop("Probability thresholds must be in [0, 1].")
    } else if (value <= 0 || value >= 1) stop("Probability mean/mode targets must be strictly inside (0, 1).")
  } else if (value < 0 || (record$type != "probability" && value <= 0)) stop("Expected-count targets must be positive (inequality thresholds may be zero).")
  record
}

prior_outcome_link = function(value, scale) {
  if (scale == "probability") stats::qlogis(value) else log(value)
}

prior_outcome_probability_records = function(records) {
  lapply(records, prior_outcome_probability_record)
}
prior_outcome_probability_record = function(record) {
  if (prior_outcome_scale(record$scale) && record$type == "probability") record$threshold = prior_outcome_link(record$threshold, record$scale)
  record
}

prior_outcome_integrand = function(z, location, scale, df, power) {
  value = stats::plogis(location + scale * z)^power
  value * if (is.null(df)) stats::dnorm(z) else stats::dt(z, df)
}

prior_logit_normal_mode = function(location, scale) {
  variance = scale^2
  bound = abs(location) + variance + 50
  points = c(-bound, bound)
  if (variance > 2) {
    critical = stats::qlogis((1 + sqrt(1 - 2 / variance)) / 2)
    points = c(-bound, -critical, critical, bound)
  }

  # Search every monotone segment so a second mode is not missed.
  roots = numeric()
  for (i in seq_len(length(points) - 1L)) {
    left = points[i]
    right = points[i + 1L]
    fl = prior_logit_mode_equation(left, location, variance)
    fr = prior_logit_mode_equation(right, location, variance)
    if (fl * fr <= 0) roots = c(roots, stats::uniroot(prior_logit_mode_equation, c(left, right), location = location, variance = variance, tol = 1e-12)$root)
  }
  roots = unique(roots)

  # log Jacobian evaluated stably, even near zero/one.
  density = stats::dnorm(roots, location, scale, log = TRUE) +
    pmax(roots, 0) + log1p(exp(-abs(roots))) + pmax(-roots, 0) + log1p(exp(-abs(roots)))
  winners = which(density >= max(density) - 1e-9)
  if (length(winners) != 1L) stop("Induced logit-normal has multiple global modes; a single mode is unsupported.")
  stats::plogis(roots[winners])
}
prior_logit_mode_equation = function(eta, location, variance) -(eta - location) / variance + 2 * stats::plogis(eta) - 1

prior_outcome_target = function(record, prior) {
  prior_validate_distribution(prior)
  if (record$type == "probability") {
    threshold = prior_outcome_link(record$threshold, record$scale)
    p = prior_distribution_operation(prior, "cdf", threshold)
    return(if (record$operator %in% c("<", "<=")) p else 1 - p)
  }
  normal = prior$distribution == "normal"
  location = if (normal) prior$mean else prior$location
  scale = if (normal) prior$sd else prior$scale
  if (record$scale != "probability") {
    if (!normal) stop("Exponentiated Student's t has no finite positive moments and no attained positive global mode; this belief is unsupported.")
    logvalue = switch(record$type, mean = location + scale^2 / 2,
      sd = location + scale^2 / 2 + .5 * prior_logexpm1(scale^2), mode = location - scale^2)
    value = exp(logvalue)
    if (!is.finite(value) || value <= 0) stop("Induced lognormal moment/mode is outside reliable numerical range.")
    return(value)
  }
  if (record$type == "mode") {
    if (!normal) stop("Logit-Student's t density diverges at both boundaries; an interior global mode is unsupported.")
    return(prior_logit_normal_mode(location, scale))
  }

  # Probability-scale moments need integration, with a separate error check.
  df = if (normal) NULL else prior$df
  first = stats::integrate(
    prior_outcome_integrand,
    -Inf,
    Inf,
    location = location,
    scale = scale,
    df = df,
    power = 1,
    rel.tol = 1e-10,
    abs.tol = 1e-11,
    subdivisions = 500L
  )
  if (first$abs.error > 1e-8) stop("Numerical integration error exceeds the outcome evaluation budget.")
  if (record$type == "mean") return(first$value)
  second = stats::integrate(
    prior_outcome_integrand,
    -Inf,
    Inf,
    location = location,
    scale = scale,
    df = df,
    power = 2,
    rel.tol = 1e-10,
    abs.tol = 1e-11,
    subdivisions = 500L
  )
  variance = second$value - first$value^2
  if (variance <= 0 ||
      (second$abs.error + 2 * abs(first$value) * first$abs.error) / (2 * sqrt(variance)) > 1e-8) {
    stop("Outcome SD numerical error is too large; unsupported at this precision.")
  }
  sqrt(variance)
}

prior_outcome_distribution = function(theta, family) {
  if (family == "normal") list(distribution = family, mean = theta[1], sd = exp(theta[2])) else list(
    distribution = family,
    location = theta[1],
    scale = exp(theta[2]),
    df = prior_student_t_policy()$df
  )
}
prior_outcome_values = function(theta, records, family) {
  prior = prior_outcome_distribution(theta, family)
  vapply(records, prior_distribution_target, numeric(1), prior = prior)
}
prior_outcome_objective = function(theta, records, family) {
  values = tryCatch(prior_outcome_values(theta, records, family), error = function(e) rep(NA_real_, length(records)))
  if (any(!is.finite(values))) return(1e100)
  residual = (values - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  sum(pmin(residual^2, 1e90))
}

prior_fit_outcome_unique = function(records, reference, family) {
  if (any(vapply(
      records,
      function(r) identical(r$scale, "probability") && (r$type %in% c("mean", "sd") && r$tolerance < 1e-7 || r$type == "mode" && r$tolerance < 1e-10),
      logical(1)
      ))) {
    return(list(
      valid = FALSE,
      message = "Requested tolerance is below supported numerical precision: probability mean/SD require >= 1e-7 and mode >= 1e-10; numerical evaluation has a separate 1e-8 error budget."
    ))
  }
  if (family == "student_t" && any(vapply(
      records,
      function(r) prior_outcome_scale(r$scale) && (r$scale != "probability" && r$type != "probability" || r$type == "mode"),
      logical(1)
      ))) {
    return(list(
      valid = FALSE,
      message = "Transformed Student's t has no attained interior global mode; exp(t) also has no finite positive moments. Requested combination is unsupported."
    ))
  }

  # Separate transformed moments from constraints that stay linear on the link.
  nonlinear = any(vapply(records, function(r) prior_outcome_scale(r$scale) && r$type != "probability", logical(1)))
  if (!nonlinear) {
    transformed = prior_outcome_probability_records(records)
    transformed = Filter(function(r) r$type != "probability" || is.finite(r$threshold), transformed)
    for (i in seq_along(transformed)) if (prior_outcome_scale(transformed[[i]]$scale)) transformed[[i]]$scale = "original"
    if (any(vapply(transformed, function(r) r$type == "probability" && r$target %in% c(0, 1), logical(1)))) {
      return(list(
        valid = FALSE,
        message = "A supported continuous link prior cannot attain probability zero/one at a finite interior threshold."
      ))
    }
    if (length(transformed)) {
      fit = prior_fit_location_scale(transformed, reference, family)
      if (!is.null(fit$distribution)) {
        fit$achieved = vapply(records, prior_distribution_target, numeric(1), prior = fit$distribution)
        fit$residuals = (fit$achieved - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
      }
      return(fit)
    }
    distribution = if (family == "normal") reference else list(
      distribution = family,
      location = reference$mean,
      scale = reference$sd / prior_distribution_registry()[[family]]$sd_factor,
      df = prior_student_t_policy()$df
    )
    return(list(valid = TRUE, distribution = distribution, reference_used = TRUE, score = 0))
  }
  base_scale = reference$sd / prior_distribution_registry()[[family]]$sd_factor
  starts = expand.grid(location = c(reference$mean, 0, 1, -1), scale = c(base_scale, .25, 1))
  if (length(records) == 1L) {
    if (prior_ratio_scale(records[[1]]$scale)) return(prior_fit_single_ratio(records[[1]], reference))

    # Resolve one constraint by anchoring the other parameter to a disclosed
    # reference; never alter the backend spread to make a mode feasible.
    if (records[[1]]$type == "mode" && records[[1]]$scale == "probability") return(prior_fit_single_probability_mode(records[[1]], reference))
    fixed = c(reference$mean, log(base_scale))
    index = if (records[[1]]$type == "sd") 2L else 1L
    single = stats::optimize(
      prior_outcome_single_objective,
      c(if (index == 1L) -50 else -10, if (index == 1L) 50 else 3),
      records = records,
      family = family,
      fixed = fixed,
      index = index,
      tol = 1e-12
    )
    fixed[index] = single$minimum
    values = tryCatch(prior_outcome_values(fixed, records, family), error = function(e) NA_real_)
    residuals = (values - records[[1]]$target) / records[[1]]$tolerance
    valid = all(is.finite(residuals)) && all(abs(residuals) <= 1)
    return(list(
      valid = valid,
      distribution = prior_outcome_distribution(fixed, family),
      achieved = values,
      residuals = residuals,
      score = sum(residuals^2),
      reference_used = TRUE,
      numerical = "Single-belief search: actual backend reference location or SD held fixed.",
      message = if (!valid) "No acceptable candidate at the reference anchoring; supply more beliefs." else NULL
    ))
  }
  best = NULL
  for (i in seq_len(nrow(starts))) {
    start = c(starts$location[i], log(starts$scale[i]))
    fit = stats::optim(start, prior_outcome_objective, records = records, family = family,
      method = "Nelder-Mead", control = list(maxit = 1500, reltol = 1e-13))
    if (is.null(best) || fit$value < best$value) best = fit
  }
  achieved = tryCatch(prior_outcome_values(best$par, records, family), error = function(e) e)
  if (inherits(achieved, "error")) return(list(valid = FALSE, message = conditionMessage(achieved)))
  jacobian = tryCatch(prior_outcome_jacobian(best$par, records, family), error = function(e) NULL)
  if (is.null(jacobian) || any(!is.finite(jacobian)) || min(svd(jacobian)$d) <= 1e-5 * max(svd(jacobian)$d)) {
    return(list(
      valid = FALSE,
      message = "Underidentified nonlinear constraints: automatic reference resolution for this combination is unsupported; add an SD or an informative quantile belief."
    ))
  }
  residuals = (achieved - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  valid = all(is.finite(residuals)) && all(abs(residuals) <= 1)
  list(
    valid = valid,
    distribution = prior_outcome_distribution(best$par, family),
    achieved = achieved,
    residuals = residuals,
    score = sum(residuals^2),
    reference_used = FALSE,
    numerical = if (any(vapply(records, function(r) identical(r$scale, "probability"), logical(1)))) "Outcome moments: estimated integration error budget 1e-8; mode root tolerance 1e-12; search convergence reported separately." else "Outcome moments: analytic lognormal/link evaluations; no numerical integration; search convergence reported separately.",
    convergence = best$convergence,
    message = if (!valid) "No acceptable outcome-scale candidate found; nonlinear search is not a proof of infeasibility." else NULL
  )
}

prior_outcome_single_objective = function(value, records, family, fixed, index) {
  fixed[index] = value
  prior_outcome_objective(fixed, records, family)
}

prior_outcome_constraint_key = function(record) {
  paste(c(record$type, record$scale, format(record$target, digits = 17),
    format(record$threshold, digits = 17), record$operator,
    format(record$tolerance, digits = 17)), collapse = "|")
}
prior_fit_outcome = function(records, reference, family) {
  keys = vapply(records, prior_outcome_constraint_key, character(1))
  result = prior_fit_outcome_unique(records[!duplicated(keys)], reference, family)
  if (!is.null(result$distribution) && isTRUE(result$valid)) {
    result$achieved = vapply(records, prior_distribution_target, numeric(1), prior = result$distribution)
    result$residuals = (result$achieved - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
    result$score = sum(result$residuals^2)
    result$valid = all(is.finite(result$residuals)) && all(abs(result$residuals) <= 1)
  }
  result
}

prior_outcome_jacobian = function(theta, records, family) {
  columns = list()
  for (i in 1:2) {
    left = right = theta
    left[i] = left[i] - 1e-4
    right[i] = right[i] + 1e-4
    columns[[i]] = (prior_outcome_values(right, records, family) - prior_outcome_values(left, records, family)) / 2e-4
  }
  do.call(cbind, columns)
}

prior_logexpm1 = function(x) {
  if (x > 50) x + log1p(-exp(-x)) else log(expm1(x))
}

# Stationarity of log f_P(p), including its Jacobian, fixes location uniquely
# once the requested probability mode and the reference SD are fixed.
prior_fit_single_probability_mode = function(record, reference) {
  target = record$target
  eta = stats::qlogis(target)
  location = eta - reference$sd^2 * (2 * target - 1)
  distribution = list(distribution = "normal", mean = location, sd = reference$sd)
  reason = "The actual backend-default spread cannot support this requested mode as a unique interior global mode (or cannot establish it reliably near a numerical boundary). Add an SD or informative quantile belief; no reference spread was modified."
  failure = list(valid = FALSE, reference_used = FALSE, reference_sd = reference$sd, message = reason)

  # Positive curvature makes this stationary point a local minimum.
  if (-1 / reference$sd^2 + 2 * target * (1 - target) > 1e-12) return(failure)
  achieved = tryCatch(prior_logit_normal_mode(location, reference$sd), error = function(e) NA_real_)
  if (!is.finite(achieved) || achieved <= 0 || achieved >= 1 ||
      abs(stats::qlogis(achieved) - eta) > 1e-8 * (1 + abs(eta))) {
    return(failure)
  }
  residual = (achieved - target) / record$tolerance
  if (abs(residual) > 1) return(failure)
  list(
    valid = TRUE,
    distribution = distribution,
    achieved = achieved,
    residuals = residual,
    score = residual^2,
    reference_used = TRUE,
    reference_sd = reference$sd,
    numerical = paste0("Single probability mode: actual backend-default link SD = ",
      format(reference$sd), " preserved unchanged; location elicited using the Jacobian and unique global mode verified. Mode root tolerance 1e-12.")
  )
}
