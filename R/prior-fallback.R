# Selection history is ordinary R data: updates never mutate an earlier object.
prior_retention_initialize = function(object) {
  if (is.null(object$fallbacks)) object$fallbacks = list()
  if (is.null(object$retention_model)) object$retention_model = prior_retention_model(object)
  if (!is.null(object$last_valid)) return(object)
  object$last_valid = object$selected
  object$last_valid_source = stats::setNames(rep("backend_default", length(object$selected)), names(object$selected))
  for (name in names(object$selected)) {
    if (is.null(object$last_valid[[name]])) object$last_valid[name] = list(object$defaults[[name]]) else
      if (length(object$beliefs[[name]]) || (!is.null(object$joint) &&
          !identical(object$selected[[name]], object$defaults[[name]]))) object$last_valid_source[name] = "elicited"
  }
  object
}

prior_retention_model = function(object) {
  list(formula = object$model$formula, family = object$model$family$family,
    link = object$model$family$link, reference = object$model$reference,
    contrasts = object$model$contrasts, package = object$backend$package,
    version = object$backend$version, QR = object$backend$QR, sparse = object$backend$sparse)
}

prior_available_compatible = function(object, selected = object$selected) {
  if (any(vapply(selected, is.null, logical(1)))) return(FALSE)
  valid = vapply(selected, function(prior) tryCatch({ prior_validate_distribution(prior); TRUE }, error = function(e) FALSE), logical(1))
  if (!all(valid)) return(FALSE)
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  length(unique(vapply(selected[slopes], `[[`, character(1), "distribution"))) <= 1L
}

prior_fallback_record = function(object, name, status, message) {
  prior = object$last_valid[[name]]
  records = object$beliefs[[name]]
  achieved = tryCatch(vapply(records, prior_distribution_target, numeric(1), prior = prior), error = function(e) rep(NA_real_, length(records)))
  within = rep(NA, length(records))
  evaluated = is.finite(achieved)
  targets = vapply(records, `[[`, numeric(1), "target")
  tolerances = vapply(records, `[[`, numeric(1), "tolerance")
  within[evaluated] = abs(achieved[evaluated] - targets[evaluated]) <= tolerances[evaluated]
  list(prior = prior, source = unname(object$last_valid_source[name]), status = status,
    reason = message, beliefs = records, verified = FALSE,
    satisfied = if (length(within) && any(!within, na.rm = TRUE)) FALSE else NA,
    achieved = achieved, within_tolerance = within)
}

prior_record_success = function(object) {
  for (name in names(object$selected)) {
    if (!is.null(object$fallbacks[[name]])) next
    object$last_valid[name] = list(object$selected[[name]])
    elicited = length(object$beliefs[[name]]) || (!is.null(object$joint) &&
      !identical(object$selected[[name]], object$defaults[[name]]))
    object$last_valid_source[name] = if (elicited) "elicited" else "backend_default"
  }
  object
}

prior_select_compatible = function(object) {
  object = prior_retention_initialize(object)
  previous = object$selected
  previous_fallbacks = object$fallbacks
  object = prior_choose_candidates(object)
  if (!length(object$outcomes)) {
    object$joint = NULL
    return(prior_record_success(object))
  }

  # An unresolved coefficient must not become a soft constraint in reconciliation.
  blocked = intersect(names(object$fallbacks), colnames(object$model$X))
  unsupported = any(vapply(object$selected[colnames(object$model$X)], function(prior) prior$distribution != "normal", logical(1)))
  if (length(blocked) || unsupported) {
    design = prior_joint_design(object, object$outcomes)
    message = if (length(blocked)) paste("Outcome reconciliation not attempted: resolve coefficient elicitation for", paste(blocked, collapse = ", "), ".") else
      "Outcome reconciliation not attempted: it requires Normal coefficient priors; Student's t coefficient priors are unsupported. Mathematical incompatibility was not established."
    object$joint = list(status = "unresolved", beliefs = object$beliefs,
      specifications = object$outcomes, baseline = object$selected,
      attempted = FALSE, design = design, message = message)
  } else object = prior_joint_reconcile(object)
  if (object$joint$status %in% c("compatible", "reconciled")) return(prior_record_success(object))

  # Never assemble a failed joint proposal from pieces of different solutions.
  if (!prior_available_compatible(object, previous)) stop("No coherent retained prior collection is available; identify priors again.")
  proposed = object$selected
  object$selected = previous
  if (prior_has_sigma(object)) object$selected$sigma = proposed$sigma
  object$joint$retained = TRUE
  object$joint$retained_selected = object$selected
  affected = colnames(object$model$X)
  for (name in affected) {
    object$last_valid[name] = list(previous[[name]])
    fallback = prior_fallback_record(object, name, object$joint$status, object$joint$message)
    if (!is.null(previous_fallbacks[[name]])) fallback$source = previous_fallbacks[[name]]$source
    object$fallbacks[[name]] = fallback
  }
  # A successful independent sigma update survives a failed outcome proposal.
  object$backend$compatible = prior_available_compatible(object)
  prior_record_success(object)
}

prior_fallback_names = function(object) names(Filter(Negate(is.null), object$fallbacks))

prior_fallback_issues = function(object) {
  names = prior_fallback_names(object)
  issues = vapply(names, function(name) paste0(name, " (", object$fallbacks[[name]]$status, ")"), character(1))
  if (isTRUE(object$joint$retained)) issues = c(issues, paste0("outcome beliefs (", object$joint$status, ")"))
  issues
}

prior_fallback_permission = function(object, allow_fallback) {
  if (!is.logical(allow_fallback) || length(allow_fallback) != 1L || is.na(allow_fallback)) stop("allow_fallback must be TRUE or FALSE.")
  issues = prior_fallback_issues(object)
  if (length(issues) && !allow_fallback) stop("Unsuccessful elicitation: ", paste(issues, collapse = ", "),
    ". Retained priors are available, but current beliefs are unsatisfied or unverified. Revise/remove beliefs, or explicitly use allow_fallback = TRUE.", call. = FALSE)
  invisible(issues)
}

prior_validate_retained = function(object) {
  if (!identical(object$retention_model, prior_retention_model(object))) stop("Model parameterization or backend changed since prior identification; identify priors again.")
  for (name in prior_fallback_names(object)) {
    retained = object$fallbacks[[name]]
    if (!identical(object$selected[[name]], retained$prior) ||
        !identical(object$beliefs[[name]], retained$beliefs)) stop("Retained prior or belief metadata changed for ", name, "; update the beliefs again.")
  }
  if (isTRUE(object$joint$retained)) {
    joint = object$joint
    if (joint$status %in% c("compatible", "reconciled") ||
        !identical(object$selected, joint$retained_selected) ||
        !identical(object$beliefs, joint$beliefs) ||
        !identical(object$outcomes, joint$specifications)) stop("Retained outcome reconciliation metadata changed; re-elicit the priors.")
  }
  # Opt-in approves known unsatisfied beliefs, never changed model data.
  prior_validate_standardization_data(object)
  invisible(object)
}

prior_warn_fallback = function(object) {
  issues = prior_fallback_issues(object)
  if (length(issues)) warning("Using retained fallback priors for ", paste(issues, collapse = ", "),
    ". Current beliefs were not satisfied or verified; no beliefs have been discarded.", call. = FALSE)
  invisible(object)
}

prior_fallback_summary = function(object, name, diagnostics) {
  retained = object$fallbacks[[name]]
  if (is.null(retained)) return(invisible(NULL))
  origin = if (retained$source == "backend_default") "The original rstanarm default has been retained." else
    "The previous successfully selected prior has been retained."
  reason = if (diagnostics || retained$status == "contradictory") retained$reason else
    switch(retained$status,
      no_supported_fit = paste("No supported prior was verified to satisfy all current beliefs for", name, "."),
      backend_incompatible = if (name == "sigma" && prior_has_sigma(object))
        "The supported sigma families have mode zero; the requested mode is backend-incompatible." else
        "rstanarm requires one shared slope-prior family; the current candidate cannot be combined with the other priors.",
      incompatible = "Outcome beliefs are mathematically incompatible with the supported model.",
      "Current elicitation could not be verified.")
  prior_console_line(paste("Warning (", name, "): ", reason, sep = ""))
  prior_console_line(origin)
  prior_console_line(if (identical(retained$satisfied, FALSE)) "The retained prior does not satisfy all current beliefs." else
    "The retained prior is not a verified solution of the current elicitation.")
  if (!diagnostics) return(invisible(NULL))
  records = object$beliefs[[name]]
  if (length(records)) print(data.frame(Belief = seq_along(records),
    Requested = vapply(records, `[[`, numeric(1), "target"), Retained = retained$achieved,
    Within_tolerance = retained$within_tolerance), row.names = FALSE, digits = 5)
  prior_console_line("Candidate diagnostics below concern attempted fits, not an approved replacement for the retained prior.")
  invisible(NULL)
}
