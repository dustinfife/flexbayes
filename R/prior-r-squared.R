# R-squared is a derived joint quantity, never a marginal property of sigma.
prior_r_squared_consistency = function(records) {
  for (scale in c("original", "r_squared")) {
    same = records
    for (i in seq_along(same)) if (!identical(same[[i]]$scale, scale)) same[[i]]$type = "ignored"
    check = if (scale == "original") prior_scale_consistency(same) else prior_check_probability_constraints(same)
    if (check$contradictory) return(check)
  }
  list(contradictory = FALSE)
}

# Finish coefficient/outcome reconciliation first, then resolve the sigma
# dependency exactly once. Failed refits retain the entire previous collection.
prior_select_compatible = function(object) {
  previous = object
  result = prior_select_base(object)
  if (!prior_has_r_squared(result)) { result$r_squared = NULL; return(result) }
  records = result$beliefs$sigma
  context = prior_r_squared_context(result)
  consistency = prior_r_squared_consistency(records)
  blocked = intersect(prior_fallback_names(result), colnames(result$model$X))
  blocked_joint = !is.null(result$joint) && !result$joint$status %in% c("compatible", "reconciled")
  candidate = NULL
  zero = !length(context$priors)
  for (i in which(vapply(records, prior_is_r_squared, logical(1)))) {
    value = prior_r_squared_constant(records[[i]], zero)
    if (!is.null(value) && abs(value - records[[i]]$target) > records[[i]]$tolerance) {
      candidate = list(valid = FALSE, status = "contradictory",
        message = paste0("Belief ", i, " is mathematically incompatible with the R-squared support: this event has probability ", value, "."))
      break
    }
  }
  if (consistency$contradictory) candidate = list(valid = FALSE, status = "contradictory", message = consistency$message)
  if (is.null(candidate) && (length(blocked) || blocked_joint)) candidate = list(valid = FALSE,
    status = "no_supported_fit", message = "R-squared elicitation requires verified current coefficient priors; resolve coefficient/outcome beliefs first.")
  cached = previous$r_squared
  if (is.null(candidate) && !is.null(cached) && cached$status == "verified" &&
      identical(records, cached$records) && identical(context, cached$context) &&
      cached$sigma$distribution %in% prior_scale_families(result$backend$package)) {
    result$selected$sigma = cached$sigma
    result$candidates$sigma = previous$candidates$sigma
    result$fallbacks$sigma = NULL
    result$r_squared = cached
    if (!is.null(result$joint)) result$joint$adjusted = result$selected
    return(prior_record_success(result))
  }
  if (is.null(candidate)) candidate = prior_fit_scale_candidates(records, result$defaults$sigma, result$backend$package, context)
  if (!isTRUE(candidate$valid) && identical(candidate$status, "no_supported_fit") &&
      !length(blocked) && !blocked_joint && result$backend$package == "rstanarm") {
    alternative = prior_fit_scale_candidates(records, result$defaults$sigma, "brms", context)
    if (isTRUE(alternative$valid)) {
      candidate$status = "backend_incompatible"
      candidate$message = "R-squared and direct sigma beliefs have a verified brms-family candidate, but no supported rstanarm sigma candidate was verified."
      candidate$alternatives = alternative$alternatives
      candidate$mathematical = alternative$distribution
    }
  }
  result$candidates$sigma = candidate
  result$r_squared = list(status = if (isTRUE(candidate$valid)) "verified" else "unresolved",
    records = records, context = context, sigma = candidate$distribution,
    verification = candidate$verification, message = candidate$message)
  if (isTRUE(candidate$valid)) {
    result$selected$sigma = candidate$distribution
    result$fallbacks$sigma = NULL
    if (!is.null(result$joint)) result$joint$adjusted = result$selected
    return(prior_record_success(result))
  }

  # A failed dependency must never combine new coefficients with an old sigma
  # and present that collection as a verified R-squared elicitation.
  result$selected = previous$selected
  result$last_valid = previous$selected
  result$last_valid_source = previous$last_valid_source
  for (name in names(result$selected)) result$fallbacks[[name]] = prior_fallback_record(result,
    name, if (name == "sigma") candidate$status else "r_squared_unresolved", candidate$message)
  if (!is.null(result$joint)) {
    if (result$joint$status %in% c("compatible", "reconciled")) {
      result$joint$status = "unresolved"
      result$joint$message = candidate$message
    }
    result$joint$retained = TRUE
    result$joint$retained_selected = result$selected
  }
  result$r_squared$retained_coefficients = result$selected[colnames(context$covariance)]
  result$r_squared$retained_sigma = result$selected$sigma
  result$backend$compatible = prior_available_compatible(result)
  if (candidate$status != "backend_incompatible" && !blocked_joint) warning(candidate$message, call. = FALSE)
  result
}

prior_r_squared_validate = function(object, allow_fallback = FALSE) {
  if (!prior_has_r_squared(object)) return(invisible(NULL))
  context = prior_r_squared_context(object)
  info = object$r_squared
  if (allow_fallback && length(prior_fallback_issues(object))) return(invisible(NULL))
  if (is.null(info) || !identical(info$status, "verified") ||
      !identical(info$records, object$beliefs$sigma) || !identical(info$context, context) ||
      !identical(info$sigma, object$selected$sigma) || !isTRUE(info$verification$valid)) {
    stop("R-squared elicitation is unresolved or stale; update the beliefs to reevaluate sigma with the current coefficient priors.")
  }
  values = info$verification$achieved
  error = info$verification$error
  if (length(values) != length(info$records) || length(error) != length(values) ||
      any(!is.finite(c(values, error))) || any(error < 0) ||
      any(abs(values - vapply(info$records, `[[`, numeric(1), "target")) + error >
        vapply(info$records, `[[`, numeric(1), "tolerance"))) stop("R-squared probability verification metadata changed; update the beliefs again.")
  invisible(NULL)
}

prior_r_squared_summary = function(object, diagnostics = FALSE) {
  if (!prior_has_r_squared(object)) return(invisible(NULL))
  info = object$r_squared
  ready = tryCatch({ prior_r_squared_validate(object); TRUE }, error = function(e) FALSE)
  records = object$beliefs$sigma
  indices = which(vapply(records, prior_is_r_squared, logical(1)))
  rows = lapply(indices, function(i) data.frame(Belief = prior_joint_summary_belief(records[[i]], "R-squared"),
    Requested = prior_joint_summary_number(records[[i]]$target, TRUE),
    Achieved = if (ready) prior_joint_summary_number(info$verification$achieved[i], TRUE) else "Not verified"))
  cat("\nBeliefs about model R-squared\n\n")
  print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
  cat("Selected sigma prior:", prior_format(object$selected$sigma), if (!ready) "(retained; not verified)", "\n")
  prior_console_line("Explained variation under the empirical predictor distribution and current coefficient priors; individual coefficient uncertainty is integrated, and the intercept contributes no variation.")
  if (!ready && !is.null(info$message)) prior_console_line(info$message)
  if (!diagnostics || is.null(info)) return(invisible(NULL))
  cat(if (ready) "\nCoefficient priors used for verified R-squared integration:\n" else "\nAttempted coefficient priors (not selected):\n")
  for (name in names(info$context$coefficients)) cat(name, ":", prior_format(info$context$coefficients[[name]]), "\n")
  cat("Empirical model-matrix covariance (sample denominator n - 1):\n")
  print(info$context$covariance)
  if (!is.null(info$verification)) {
    cat(if (ready) "Probability verification (estimated quadrature error plus bounded omitted tail):\n" else "Attempted probabilities (not verified for retained priors):\n")
    print(data.frame(Belief = seq_along(records), Achieved = info$verification$achieved,
      Numerical_error = info$verification$error), row.names = FALSE)
    cat("Verification method:", info$verification$method, "\n")
    cat("Quadrature order:", info$verification$order, "; tail bound:", info$verification$tail_bound,
      "; mass discrepancy:", info$verification$mass_error, "\n")
  }
  cat("Numerical policy: spherical deterministic quadrature, common nodes, analytic conditional sigma probabilities; finite order/point/dimension budgets. Error estimates are not formal global error certificates.\n")
  candidate = object$candidates$sigma
  if (!is.null(candidate$numerical)) prior_console_line(candidate$numerical)
  for (family in names(candidate$alternatives)) {
    value = candidate$alternatives[[family]]
    prior_console_line(paste(prior_candidate_label(family), if (isTRUE(value$valid)) "verified" else value$message))
  }
  # This observed-outcome comparison is diagnostic only, never a fitting input.
  means = vapply(info$context$coefficients, prior_distribution_operation, numeric(1), operation = "mean")
  variances = vapply(info$context$coefficients, prior_distribution_operation, numeric(1), operation = "sd")^2
  S = info$context$covariance
  explained = sum(S * outer(means, means)) + sum(diag(S) * variances)
  cat(if (ready) "Mean implied explained variance:" else "Attempted mean explained variance (not selected):", prior_scale_number(explained),
    "; observed outcome variance:", prior_scale_number(info$context$observed_variance), "\n")
  if (is.finite(explained) && explained > info$context$observed_variance) prior_console_line("Diagnostic discrepancy: mean implied explained variance exceeds observed outcome variance. Observed outcome variance was not added as a constraint.")
  for (record in Filter(prior_is_r_squared, records)) {
    h = record$threshold
    if (h <= 0 || h >= 1) next
    prior_console_line(paste0("At the belief threshold R-squared = ", format(h),
      ", the observed-variance shortcut suggests sigma = ", prior_scale_number(sqrt(info$context$observed_variance * (1 - h))),
      " and explained variance = ", prior_scale_number(info$context$observed_variance * h),
      ". These are diagnostic benchmarks at a threshold, not means or additional constraints."))
  }
  prior_console_line("This observed-variance comparison is not a mathematical incompatibility certificate.")
  prior_console_line("SD(Y)*sqrt(1-R-squared) is an observed-variance shortcut, not the definition or a sigma-fitting constraint used here.")
  invisible(NULL)
}

prior_r_squared_plot = function(object) {
  value = probability = NULL
  prior_predictive_validate(object)
  context = prior_r_squared_context(object)
  order = if (is.null(object$r_squared)) 64L else max(64L, object$r_squared$verification$order)
  grid = prior_r_squared_grid(context, order)
  values = seq(0, 1, length.out = 101L)
  cdf = vapply(values, function(h) prior_r_squared_probability(list(threshold = h,
    operator = "<=", type = "probability"), object$selected$sigma, grid), numeric(1))
  ggplot2::ggplot(data.frame(value = values, probability = pmin(1, pmax(0, cdf))), ggplot2::aes(x = value, y = probability)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::labs(x = "Model R-squared",
      y = "Cumulative prior probability", title = "Prior distribution of model R-squared",
      caption = "Numerical CDF under the recorded empirical predictor distribution.")
}
