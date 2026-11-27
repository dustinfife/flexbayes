prior_r_squared_cached_grid = function(context, order) {
  key = as.character(order)
  if (!is.null(context$cache) && exists(key, context$cache, inherits = FALSE)) return(get(key, context$cache))
  grid = prior_r_squared_grid(context, order)
  if (!is.null(context$cache)) grid$budget = context$cache
  if (!is.null(context$cache)) assign(key, grid, context$cache)
  grid
}

prior_r_squared_distribution = function(family, scale, shape) {
  if (family %in% c("gamma", "lognormal")) return(prior_positive_distribution(family, shape, scale))
  prior_scale_distribution(family, scale, if (family == "half_student_t") shape else NULL)
}

prior_r_squared_scale_residual = function(log_scale, family, shape, record, grid, reference) {
  prior = prior_r_squared_distribution(family, exp(log_scale) / reference$rate, shape)
  value = if (prior_is_r_squared(record)) prior_r_squared_probability(record, prior, grid) else prior_distribution_target(record, prior)
  value - record$target
}

prior_r_squared_profile = function(log_shape, family, records, reference, grid, floor = 0) {
  shape = if (family %in% c("half_student_t", "gamma", "lognormal")) floor + exp(log_shape) else NULL
  informative = which(vapply(records, function(r) prior_is_r_squared(r) && r$threshold > 0 && r$threshold < 1 ||
    !prior_is_r_squared(r) && (r$type %in% c("mean", "sd") || r$type == "mode" && r$target > 0 || r$type == "probability" && r$threshold > 0), logical(1)))
  if (!length(informative)) return(list(distribution = prior_r_squared_distribution(family, 1 / reference$rate, shape), residual = rep(0, length(records))))
  direct = informative[!vapply(records[informative], prior_is_r_squared, logical(1))]
  anchor = prior_r_squared_anchor(records[[if (length(direct)) direct[1] else informative[1]]], records)
  unit = prior_r_squared_distribution(family, 1, shape)
  if (!prior_is_r_squared(anchor) && anchor$type != "probability") {
    scale = anchor$target / prior_distribution_operation(unit, anchor$type)
  } else {
    # Fixed coefficient nodes are shared by every scale and shape trial.
    root = stats::uniroot(prior_r_squared_scale_residual, c(-30, 30),
      family = family, shape = shape, record = anchor, grid = grid, reference = reference, tol = 1e-11)
    scale = exp(root$root) / reference$rate
  }
  prior = prior_r_squared_distribution(family, scale, shape)
  prior_validate_distribution(prior)
  achieved = prior_r_squared_values(prior, records, grid)
  units = vapply(records, function(r) if (r$type == "probability") 1 else 1 / reference$rate, numeric(1))
  list(distribution = prior, residual = (achieved - vapply(records, `[[`, numeric(1), "target")) / units)
}

prior_r_squared_profile_score = function(log_shape, ...) {
  result = tryCatch(prior_r_squared_profile(log_shape, ...), error = function(e) NULL)
  if (is.null(result) || any(!is.finite(result$residual))) return(1e100)
  sum(result$residual^2)
}
prior_r_squared_profile_residual = function(log_shape, index, ...) prior_r_squared_profile(log_shape, ...)$residual[index]

prior_fit_r_squared_family = function(family, records, reference, context) {
  direct = Filter(function(r) !prior_is_r_squared(r), records)
  positive_mode = any(vapply(direct, function(r) r$type == "mode" && r$target > r$tolerance, logical(1)))
  if (positive_mode && !family %in% c("gamma", "lognormal")) return(list(valid = FALSE, message = "This sigma family has mode zero."))
  shape_family = family %in% c("half_student_t", "gamma", "lognormal")
  informative = Filter(function(r) prior_is_r_squared(r) && r$threshold > 0 && r$threshold < 1 ||
    !prior_is_r_squared(r) && (r$type != "mode" || r$target > 0), records)
  reference_used = length(unique(vapply(informative, prior_r_squared_information_key, character(1)))) < if (shape_family) 2L else 1L
  fit_records = records
  if (reference_used && family %in% c("gamma", "lognormal")) fit_records = c(records,
    list(list(type = "sd", target = 1 / reference$rate, tolerance = 1e-6, scale = "original")))
  floor = if (family == "half_student_t") {
    if (any(vapply(direct, function(r) r$type == "sd", logical(1)))) 2 else if (any(vapply(direct, function(r) r$type == "mean", logical(1)))) 1 else 0
  } else if (family == "gamma" && positive_mode) 1 else 0
  bounds = if (family == "half_student_t") c(-14, log(1e6 - floor)) else c(-12, 12)
  reference_shape = if (family == "half_student_t") log(prior_scale_df_policy()$reference - floor) else 0
  best = NULL
  for (resolution in c(16L, 32L, 64L, 128L)) {
    grid = prior_r_squared_cached_grid(context, resolution)
    trials = if (shape_family && !(reference_used && family == "half_student_t"))
      sort(unique(c(seq(bounds[1], bounds[2], length.out = 17L), reference_shape))) else reference_shape
    arguments = list(family = family, records = fit_records, reference = reference, grid = grid, floor = floor)
    if (length(trials) > 1L) {
      scores = vapply(trials, prior_r_squared_profile_score, numeric(1), family = family,
        records = fit_records, reference = reference, grid = grid, floor = floor)
      roots = numeric()
      for (i in seq_len(length(trials) - 1L)) {
        a = tryCatch(do.call(prior_r_squared_profile, c(list(trials[i]), arguments))$residual, error = function(e) NULL)
        b = tryCatch(do.call(prior_r_squared_profile, c(list(trials[i + 1L]), arguments))$residual, error = function(e) NULL)
        if (is.null(a) || is.null(b)) next
        for (index in which(a * b < 0)) {
          root = tryCatch(do.call(stats::uniroot, c(list(f = prior_r_squared_profile_residual,
            interval = trials[c(i, i + 1L)], index = index, tol = 1e-10), arguments))$root, error = function(e) NULL)
          roots = c(roots, root)
        }
      }
      minima = which(scores <= c(Inf, head(scores, -1)) & scores <= c(tail(scores, -1), Inf))
      for (i in head(minima, 3L)) {
        fit = tryCatch(do.call(stats::optimize, c(list(f = prior_r_squared_profile_score,
          interval = trials[c(max(1L, i - 1L), min(length(trials), i + 1L))], tol = 1e-10), arguments)), error = function(e) NULL)
        if (!is.null(fit)) roots = c(roots, fit$minimum)
      }
      trials = unique(c(roots, trials[order(scores)]))
    }
    for (trial in trials) {
      fit = tryCatch(do.call(prior_r_squared_profile, c(list(trial), arguments)), error = function(e) NULL)
      if (is.null(fit)) next
      approximate = prior_r_squared_values(fit$distribution, records, grid)
      tolerance = vapply(records, `[[`, numeric(1), "tolerance")
      target = vapply(records, `[[`, numeric(1), "target")
      if (any(abs(approximate - target) > tolerance)) next
      verification = prior_r_squared_verify(fit$distribution, records, context)
      best = list(valid = verification$valid, distribution = fit$distribution,
        achieved = verification$achieved, residuals = (verification$achieved - target) / tolerance,
        score = sum(((verification$achieved - target) / tolerance)^2),
        reference_used = reference_used && family != "half_student_t",
        verification = verification, numerical = paste("Deterministic spherical Gauss-Legendre/periodic quadrature; sigma integrated analytically.",
          "Success requires discrepancy plus estimated numerical error within every original tolerance.",
          if (reference_used && family == "half_student_t") "Underdetermined half-t uses reference df = 7." else
          if (reference_used) "An underdetermined positive family references the default sigma SD." else "No reference constraint added."))
      if (isTRUE(best$valid)) return(best)
    }
  }
  if (is.null(best)) best = list(valid = FALSE)
  best$message = "No sigma candidate verified with numerical error inside the original tolerances; bounded numerical search does not prove incompatibility."
  best
}

prior_fit_r_squared_candidates = function(records, reference, backend, context) {
  r2 = Filter(prior_is_r_squared, records)
  constant = vapply(r2, function(r) !is.null(prior_r_squared_constant(r, !length(context$priors))), logical(1))
  if (all(constant)) return(prior_fit_r_squared_constants(records, reference, backend, context))
  uncertain = r2[!constant]
  if (any(vapply(uncertain, `[[`, numeric(1), "tolerance") < 8 * prior_r_squared_policy()$tail)) return(list(
    valid = FALSE, status = "no_supported_fit", alternatives = list(),
    message = "R-squared probability cannot be verified at the requested tolerance with the current numerical tail budget. No tolerance was widened."))
  context$cache = new.env(parent = emptyenv())
  context$cache$points = 0
  alternatives = list()
  found = FALSE
  for (family in prior_scale_families(backend)) {
    candidate = if (found) list(valid = FALSE, message = "Not searched: a higher-ranked family was verified.") else
      tryCatch(prior_fit_r_squared_family(family, records, reference, context),
        error = function(e) list(valid = FALSE, message = paste("Numerical R-squared search unresolved:", conditionMessage(e))))
    candidate$family = family
    alternatives[[family]] = candidate
    found = found || isTRUE(candidate$valid)
  }
  prior_scale_rank_candidates(alternatives)
}

prior_r_squared_information_key = function(record) {
  paste(record$scale, if (record$type == "probability") paste("cdf", sprintf("%.17g", record$threshold)) else record$type)
}

# Endpoint/support statements need no quadrature, even in high dimensions.
prior_fit_r_squared_constants = function(records, reference, backend, context) {
  direct = Filter(function(r) !prior_is_r_squared(r), records)
  fit = if (length(direct)) prior_fit_scale_candidates(direct, reference, backend) else
    list(valid = TRUE, status = "selected", distribution = reference, family = "exponential", reference_used = TRUE, alternatives = list())
  if (!isTRUE(fit$valid)) return(fit)
  grid = list(v = if (length(context$priors)) 1 else 0, w = 1, zero = !length(context$priors))
  values = prior_r_squared_values(fit$distribution, records, grid)
  fit$verification = list(valid = TRUE, precise = TRUE, achieved = values, error = rep(0, length(values)),
    order = 0L, method = "Deterministic support plus direct sigma constraints", mass_error = 0, tail_bound = 0)
  fit$achieved = values
  fit$residuals = (values - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  fit$valid = fit$verification$valid = all(abs(fit$residuals) <= 1)
  fit$score = sum(fit$residuals^2)
  fit
}

# Repeated statements about the same property use the intersection of their
# declared bands for initialization. Every original answer is still verified.
prior_r_squared_anchor = function(anchor, records) {
  same = Filter(function(r) identical(r$scale, anchor$scale) && r$type == anchor$type &&
    (r$type != "probability" || r$threshold == anchor$threshold), records)
  if (anchor$type != "probability") {
    lower = max(vapply(same, function(r) r$target - r$tolerance, numeric(1)))
    upper = min(vapply(same, function(r) r$target + r$tolerance, numeric(1)))
    if (lower <= upper) anchor$target = lower + (upper - lower) / 2
    return(anchor)
  }
  probabilities = vapply(same, function(r) if (r$operator %in% c("<", "<=")) r$target else 1 - r$target, numeric(1))
  tolerances = vapply(same, `[[`, numeric(1), "tolerance")
  lower = max(0, probabilities - tolerances)
  upper = min(1, probabilities + tolerances)
  if (lower > upper) return(anchor)
  probability = lower + (upper - lower) / 2
  anchor$target = if (anchor$operator %in% c("<", "<=")) probability else 1 - probability
  anchor
}
