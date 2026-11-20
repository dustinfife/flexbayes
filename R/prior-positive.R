# Positive sigma families use their actual moments, not latent-log moments.
prior_gamma_survival = function(q, shape, rate) stats::pgamma(q, shape, rate, lower.tail = FALSE)
prior_gamma_mean = function(prior) prior$shape / prior$rate
prior_gamma_sd = function(prior) sqrt(prior$shape) / prior$rate
prior_gamma_mode = function(prior) max(0, prior$shape - 1) / prior$rate
prior_lnorm_survival = function(q, meanlog, sdlog) stats::plnorm(q, meanlog, sdlog, lower.tail = FALSE)
prior_lnorm_mean = function(prior) exp(prior$meanlog + prior$sdlog^2 / 2)
prior_lnorm_sd = function(prior) exp(prior$meanlog + prior$sdlog^2 / 2) * sqrt(expm1(prior$sdlog^2))
prior_lnorm_mode = function(prior) exp(prior$meanlog - prior$sdlog^2)

prior_positive_distribution = function(family, shape, scale = 1) {
  if (family == "gamma") return(list(distribution = family, shape = shape, rate = 1 / scale))
  list(distribution = family, meanlog = log(scale), sdlog = shape)
}

prior_positive_profile = function(log_shape, family, anchor, reference, records) {
  shape = if (family == "gamma" && anchor$type == "mode" && anchor$target > 0) 1 + exp(log_shape) else exp(log_shape)
  unit = prior_positive_distribution(family, shape)
  if (anchor$type == "probability") {
    probability = anchor$target
    lower = anchor$operator %in% c("<", "<=")
    quantile = if (family == "gamma") stats::qgamma(probability, shape, lower.tail = lower) else stats::qlnorm(probability, 0, shape, lower.tail = lower)
    scale = anchor$threshold / quantile
  } else scale = anchor$target / prior_distribution_operation(unit, anchor$type)
  prior = prior_positive_distribution(family, shape, scale)
  prior_validate_distribution(prior)
  achieved = vapply(records, prior_distribution_target, numeric(1), prior = prior)
  units = vapply(records, function(r) if (r$type == "probability") 1 else 1 / reference$rate, numeric(1))
  list(distribution = prior, residual = (achieved - vapply(records, `[[`, numeric(1), "target")) / units)
}

prior_positive_score = function(log_shape, ...) {
  result = tryCatch(prior_positive_profile(log_shape, ...)$residual, error = function(e) NA_real_)
  if (any(!is.finite(result))) return(1e100)
  sum(result^2)
}

prior_positive_residual = function(log_shape, index, ...) prior_positive_profile(log_shape, ...)$residual[index]

prior_positive_constraint_key = function(record) {
  if (record$type != "probability") return(paste(record$type, sprintf("%.17g", record$target)))
  probability = if (record$operator %in% c("<", "<=")) record$target else 1 - record$target
  paste("cdf", sprintf("%.17g", record$threshold), sprintf("%.17g", probability))
}

prior_fit_positive = function(family, records, reference) {
  informative = which(vapply(records, function(r) (r$type %in% c("mean", "sd", "mode") && r$target > 0) ||
    (r$type == "probability" && r$threshold > 0 && r$target > 0 && r$target < 1), logical(1)))
  if (!length(informative)) return(list(valid = FALSE, message = "No informative positive-family anchor; exponential reference preferred."))
  anchor = records[[informative[1]]]
  distinct = unique(vapply(records[informative], prior_positive_constraint_key, character(1)))
  reference_used = length(distinct) == 1L
  fit_records = records
  if (reference_used) {
    # Preserve an actual reference moment, rather than inserting a shape or cap.
    moment = if (anchor$type == "sd") "mean" else "sd"
    fit_records = c(records, list(list(type = moment, target = 1 / reference$rate, tolerance = 1e-6, scale = "original")))
  }
  grid = seq(-12, 12, length.out = 121)
  trials = grid
  scores = vapply(grid, prior_positive_score, numeric(1), family = family, anchor = anchor, reference = reference, records = fit_records)
  for (i in seq_len(length(grid) - 1L)) {
    a = tryCatch(prior_positive_profile(grid[i], family, anchor, reference, fit_records)$residual, error = function(e) NULL)
    b = tryCatch(prior_positive_profile(grid[i + 1L], family, anchor, reference, fit_records)$residual, error = function(e) NULL)
    if (is.null(a) || is.null(b)) next
    for (index in which(is.finite(a) & is.finite(b) & a * b < 0)) {
      root = tryCatch(stats::uniroot(prior_positive_residual, grid[c(i, i + 1L)], index = index,
        family = family, anchor = anchor, reference = reference, records = fit_records, tol = 1e-12)$root, error = function(e) NULL)
      trials = c(trials, root)
    }
  }
  minima = which(scores <= c(Inf, head(scores, -1)) & scores <= c(tail(scores, -1), Inf))
  for (i in minima) {
    fit = stats::optimize(prior_positive_score, grid[c(max(1L, i - 1L), min(length(grid), i + 1L))],
      family = family, anchor = anchor, reference = reference, records = fit_records, tol = 1e-12)
    trials = c(trials, fit$minimum)
  }
  candidates = lapply(trials, function(value) {
    prior = tryCatch(prior_positive_profile(value, family, anchor, reference, fit_records)$distribution, error = function(e) NULL)
    result = prior_scale_evaluate(prior, fit_records, reference_used)
    if (isTRUE(result$valid)) result = prior_scale_evaluate(prior, records, reference_used)
    result
  })
  valid = which(vapply(candidates, function(c) isTRUE(c$valid), logical(1)))
  if (!length(valid)) return(list(valid = FALSE, message = "No positive-family candidate verified; finite log-shape search [-12,12] does not prove mathematical incompatibility."))
  best = candidates[[valid[which.min(vapply(candidates[valid], `[[`, numeric(1), "score"))]]]
  best$numerical = paste("Positive-family scale profiled from a belief; finite shape search, with all original tolerances checked.",
    if (reference_used) paste("Underdetermined: retained the universal exponential default's", moment, "as a reference.") else "No reference moment added.")
  best
}
