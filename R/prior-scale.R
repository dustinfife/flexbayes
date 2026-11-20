prior_scale_families = function(backend = "rstanarm") c("exponential", "half_normal", "half_cauchy", "half_student_t",
  if (backend == "brms") c("gamma", "lognormal"))

# This policy is separate from the fixed df = 7 policy for coefficient priors.
# Finite search limits are numerical limits, not a mathematical feasibility test.
prior_scale_df_policy = function() list(reference = 7, log_offset_min = -14, maximum = 1e6, grid = 81L)

prior_scale_consistency = function(records) {
  # Sigma is strictly positive; events below or at zero have probability zero.
  checks = records
  for (i in seq_along(checks)) {
    if (checks[[i]]$type == "probability" && checks[[i]]$threshold <= 0) checks[[i]]$threshold = -Inf
  }
  prior_check_probability_constraints(checks)
}

prior_scale_distribution = function(family, scale, df = NULL) {
  if (family == "exponential") return(list(distribution = family, rate = 1 / scale))
  if (family == "half_student_t") return(list(distribution = family, scale = scale, df = df))
  list(distribution = family, scale = scale)
}

prior_scale_anchor = function(records) {
  informative = which(vapply(records, function(r) r$type %in% c("mean", "sd") ||
    (r$type == "probability" && r$threshold > 0), logical(1)))
  if (!length(informative)) return(NULL)
  records[[informative[1]]]
}

prior_scale_from_belief = function(record, family, df, reference) {
  if (is.null(record)) return(1 / reference$rate)
  unit = prior_scale_distribution(family, 1, df)
  if (record$type != "probability") {
    value = prior_distribution_operation(unit, record$type)
    if (!is.finite(value) || value <= 0) stop("Requested moment does not exist for this family and degrees of freedom.")
    return(record$target / value)
  }

  # Boundary targets use an interior starting value of their declared band.
  # Every candidate is still checked against the original target and tolerance.
  p = record$target
  if (p == 0 || p == 1) p = (max(0, p - record$tolerance) + min(1, p + record$tolerance)) / 2
  lower = record$operator %in% c("<", "<=")
  q = prior_scale_unit_quantile(p, family, df, lower)
  record$threshold / q
}

prior_scale_unit_quantile = function(p, family, df, lower) {
  switch(family,
    exponential = stats::qexp(p, lower.tail = lower),
    half_normal = if (lower) prior_qhn(p, 1) else stats::qnorm(p / 2, lower.tail = FALSE),
    half_student_t = if (lower) prior_qht(p, 1, df) else stats::qt(p / 2, df, lower.tail = FALSE),
    half_cauchy = if (lower) prior_qhc(p, 1) else stats::qt(p / 2, 1, lower.tail = FALSE))
}

# For a fixed shape, every supported constraint gives an interval for scale.
# Intersect the actual acceptance bands, not just their point targets.
prior_fit_scale_fixed = function(family, records, reference, df = NULL) {
  unit = prior_scale_distribution(family, 1, df)
  lower = 0
  upper = Inf
  for (record in records) {
    if (record$type == "mode") next
    if (record$type == "probability") {
      if (record$threshold <= 0) next
      p = c(max(0, record$target - record$tolerance), min(1, record$target + record$tolerance))
      lower_tail = record$operator %in% c("<", "<=")
      q = prior_scale_unit_quantile(p, family, df, lower_tail)
      interval = sort(record$threshold / q)
    } else {
      moment = prior_distribution_operation(unit, record$type)
      if (!is.finite(moment) || moment <= 0) return(list(valid = FALSE,
        message = "Requested mean or SD has no finite value for this family."))
      interval = c(max(0, record$target - record$tolerance), record$target + record$tolerance) / moment
    }
    lower = max(lower, interval[1])
    upper = min(upper, interval[2])
  }
  if (lower > upper || upper <= 0 || !is.finite(lower)) return(list(valid = FALSE,
    message = "No positive scale in the intersection of this family's acceptance bands."))

  anchor = prior_scale_anchor(records)
  preferred = tryCatch(prior_scale_from_belief(anchor, family, df, reference), error = function(e) 1 / reference$rate)
  clamped = min(upper, max(lower, preferred))
  interior = if (is.finite(upper)) lower + (upper - lower) / 2 else max(1 / reference$rate, lower * 2)
  scales = unique(c(preferred, clamped, interior))
  scales = scales[is.finite(scales) & scales > 0 & scales >= lower & scales <= upper]
  candidates = lapply(scales, function(scale) prior_scale_evaluate(
    prior_scale_distribution(family, scale, df), records, is.null(anchor)))
  valid = which(vapply(candidates, function(c) isTRUE(c$valid), logical(1)))
  if (!length(valid)) return(list(valid = FALSE, message = "No candidate verified at numerical precision within the original scale bands."))
  candidates[[valid[which.min(vapply(candidates[valid], `[[`, numeric(1), "score"))]]]
}

prior_scale_profile = function(log_offset, records, reference, anchor, floor) {
  df = floor + exp(log_offset)
  scale = prior_scale_from_belief(anchor, "half_student_t", df, reference)
  prior = prior_scale_distribution("half_student_t", scale, df)
  prior_validate_distribution(prior)
  achieved = vapply(records, prior_distribution_target, numeric(1), prior = prior)
  targets = vapply(records, `[[`, numeric(1), "target")
  units = vapply(records, function(r) if (r$type == "probability") 1 else 1 / reference$rate, numeric(1))
  list(distribution = prior, residual = (achieved - targets) / units)
}

prior_scale_profile_score = function(log_offset, ...) {
  result = tryCatch(prior_scale_profile(log_offset, ...)$residual, error = function(e) NA_real_)
  if (any(!is.finite(result))) return(1e100)
  sum(result^2)
}

prior_scale_profile_residual = function(log_offset, index, ...) prior_scale_profile(log_offset, ...)$residual[index]

prior_fit_half_t = function(records, reference, anchor) {
  policy = prior_scale_df_policy()
  types = vapply(records, `[[`, character(1), "type")
  floor = if ("sd" %in% types) 2 else if ("mean" %in% types) 1 else 0
  grid = sort(unique(c(seq(policy$log_offset_min, log(policy$maximum - floor), length.out = policy$grid),
    log(policy$reference - floor))))
  scores = vapply(grid, prior_scale_profile_score, numeric(1), records = records,
    reference = reference, anchor = anchor, floor = floor)
  trials = grid

  # Profile out scale analytically, then search df without treating convergence
  # as proof that all of the original constraints are satisfied.
  for (i in seq_len(length(grid) - 1L)) {
    a = tryCatch(prior_scale_profile(grid[i], records, reference, anchor, floor)$residual, error = function(e) NULL)
    b = tryCatch(prior_scale_profile(grid[i + 1L], records, reference, anchor, floor)$residual, error = function(e) NULL)
    if (is.null(a) || is.null(b)) next
    crossing = which(is.finite(a) & is.finite(b) & a * b < 0)
    for (index in crossing) {
      root = tryCatch(stats::uniroot(prior_scale_profile_residual, grid[c(i, i + 1L)],
        index = index, records = records, reference = reference, anchor = anchor,
        floor = floor, tol = 1e-12)$root, error = function(e) NULL)
      trials = c(trials, root)
    }
  }
  minima = which(scores <= c(Inf, head(scores, -1)) & scores <= c(tail(scores, -1), Inf))
  for (i in minima) {
    interval = grid[c(max(1L, i - 1L), min(length(grid), i + 1L))]
    fit = tryCatch(stats::optimize(prior_scale_profile_score, interval, records = records,
      reference = reference, anchor = anchor, floor = floor, tol = 1e-12), error = function(e) NULL)
    if (!is.null(fit)) trials = c(trials, fit$minimum)
  }
  candidates = lapply(trials, function(value) {
    prior_fit_scale_fixed("half_student_t", records, reference, floor + exp(value))
  })
  valid = which(vapply(candidates, function(c) isTRUE(c$valid), logical(1)))
  if (!length(valid)) return(list(valid = FALSE, message = paste(
    "No half-t candidate verified within the original tolerances.",
    "Finite df search (offset >= exp(-14), df <= 1e6) is not a proof of infeasibility.")))
  best = candidates[[valid[which.min(vapply(candidates[valid], `[[`, numeric(1), "score"))]]]
  best$numerical = "Half-t scale profiled analytically; df fitted with deterministic grid/root/minimum searches. No uniqueness or global infeasibility claim. Original belief tolerances retained."
  best
}

prior_scale_evaluate = function(prior, records, reference_used = FALSE) {
  if (is.null(prior)) return(list(valid = FALSE, message = "Numerical fitting did not produce a distribution."))
  result = tryCatch({
    prior_validate_distribution(prior)
    vapply(records, prior_distribution_target, numeric(1), prior = prior)
  }, error = identity)
  if (inherits(result, "error")) return(list(valid = FALSE, message = conditionMessage(result)))
  if (any(!is.finite(result))) return(list(valid = FALSE, message = "Requested mean or SD has no finite value for this family."))
  residual = (result - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  valid = all(is.finite(residual)) && all(abs(residual) <= 1)
  list(valid = valid, distribution = prior, achieved = result, residuals = residual,
    score = sum(residual^2), reference_used = reference_used,
    message = if (!valid) "This family does not satisfy all original belief tolerances." else NULL)
}

prior_fit_scale_family = function(family, records, reference) {
  anchor = prior_scale_anchor(records)
  if (family == "half_student_t" && !is.null(anchor)) return(prior_fit_half_t(records, reference, anchor))
  result = prior_fit_scale_fixed(family, records, reference,
    if (family == "half_student_t") prior_scale_df_policy()$reference else NULL)
  result$numerical = if (is.null(anchor)) "Backend-default scale retained to resolve uninformative constraints." else
    "Scale chosen from the joint intersection of all original acceptance bands; every constraint verified. No uniqueness across families is claimed."
  result
}

prior_fit_scale_candidates = function(records, reference, backend = "rstanarm") {
  families = prior_scale_families(backend)
  incompatible = any(vapply(records, function(r) r$type == "mode" && abs(r$target) > r$tolerance, logical(1)))
  if (incompatible && backend == "rstanarm") return(list(valid = FALSE, status = "backend_incompatible",
    message = "The supported zero-centered rstanarm sigma families all have mode zero. The requested mode is outside its tolerance of zero; it was not reinterpreted as a mean.",
    alternatives = list(), mathematical = NULL))
  alternatives = stats::setNames(lapply(families, function(family) {
    if (incompatible && !family %in% c("gamma", "lognormal")) return(list(valid = FALSE, family = family, message = "This family has mode zero."))
    result = tryCatch(if (family %in% c("gamma", "lognormal")) prior_fit_positive(family, records, reference) else prior_fit_scale_family(family, records, reference),
      error = function(e) list(valid = FALSE, message = conditionMessage(e)))
    result$family = family
    result
  }), families)
  valid = which(vapply(alternatives, function(c) isTRUE(c$valid), logical(1)))

  # Prefer the default exponential, then one-parameter alternatives, then half-t.
  # This is deterministic selection, not evidence that the prior is unique.
  for (i in seq_along(valid)) alternatives[[valid[i]]]$rank = i
  best = if (length(valid)) alternatives[[valid[1]]] else list(valid = FALSE,
    message = paste("No supported sigma prior was verified within all original tolerances.",
      "A failed numerical search does not establish mathematical incompatibility.",
      paste(unique(vapply(alternatives, `[[`, character(1), "message")), collapse = " ")))
  best$status = if (length(valid)) "selected" else "no_supported_fit"
  best$alternatives = alternatives
  best$mathematical = if (length(valid)) best$distribution else NULL
  best
}

prior_scale_summary = function(object, diagnostics) {
  prior_console_line("These beliefs concern sigma itself: residual standard deviation, not residual errors or expected-response uncertainty.")
  prior = object$selected$sigma
  if (is.null(prior)) return(invisible(NULL))
  mean = prior_distribution_operation(prior, "mean")
  sd = prior_distribution_operation(prior, "sd")
  prior_console_line(paste0("Sigma prior mean: ", prior_scale_number(mean),
    "; SD: ", prior_scale_number(sd), "; median: ",
    prior_scale_number(prior_distribution_operation(prior, "quantile", .5)), "; mode: ",
    prior_scale_number(prior_distribution_operation(prior, "mode")), "."))
  if (!diagnostics) return(invisible(NULL))
  prior_console_line(paste(if (object$backend$package == "brms") "Universal flexbayes default (rstanarm-derived):" else "Original backend default:", prior_format(object$defaults$sigma)))
  if (object$backend$package == "brms") prior_console_line("After exponential and the half families, Gamma then lognormal candidates are considered. A single informative constraint uses the universal default SD as reference (or its mean for an SD belief); see the selected candidate diagnostics.")
  prior_console_line("Selection prefers the backend-default family, then fewer parameters; acceptable families need not be unique. Half-t df fitting is separate from the fixed df = 7 coefficient policy.")
  candidate = object$candidates$sigma
  if (!is.null(candidate$numerical)) prior_console_line(candidate$numerical)
  invisible(NULL)
}

prior_scale_number = function(value) {
  if (is.na(value)) return("no finite moment")
  if (!is.finite(value)) return("outside numerical range")
  format(signif(value, 5), trim = TRUE, big.mark = ",")
}

prior_candidate_label = function(family) {
  if (is.null(family)) return("Distribution")
  switch(family, normal = "Normal", student_t = "Student's t", exponential = "Exponential",
    half_normal = "Half-Normal", half_student_t = "Half-Student's t", half_cauchy = "Half-Cauchy", gamma = "Gamma", lognormal = "Lognormal", "Distribution")
}
