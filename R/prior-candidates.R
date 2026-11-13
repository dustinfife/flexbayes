prior_fit_candidates = function(records, reference) {
  families = c("normal", "student_t")
  boundary = any(vapply(records, function(r) r$type == "probability" && r$target %in% c(0, 1), logical(1)))
  alternatives = stats::setNames(lapply(families, function(family) {
    outcome = any(vapply(records, function(r) prior_outcome_scale(r$scale), logical(1)))
    candidate = if (outcome) prior_fit_outcome(records, reference, family) else if (boundary) list(valid = FALSE, message = "Boundary probability is consistent, but no supported continuous family attains it at a finite threshold.") else prior_fit_location_scale(records, reference, family)
    candidate$family = family
    candidate
  }), families)
  valid = which(vapply(alternatives, function(c) isTRUE(c$valid), logical(1)))
  # Once tolerances pass, prefer the simpler normal model. Fixed-df t has
  # equal fitted dimension but adds a tail-shape assumption; numerical noise
  # in near-zero fit scores must not decide the family.
  for (i in seq_along(valid)) alternatives[[valid[i]]]$rank = i
  best = if (length(valid)) alternatives[[valid[1]]] else list(valid = FALSE, message = paste("No supported distribution fits all beliefs within tolerance.", paste(unique(vapply(alternatives, `[[`, character(1), "message")), collapse = " ")))
  best$status = if (length(valid)) "selected" else "no_supported_fit"
  best$alternatives = alternatives
  best$mathematical = if (length(valid)) alternatives[[valid[1]]]$distribution else NULL
  best
}

# Recompute automatic selections without refitting or changing candidate ranks.
# Defaults are fixed; elicited candidates are automatic choices, not user locks.
prior_select_compatible = function(object) {
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  options = lapply(names(object$selected), function(name) {
    if (!length(object$beliefs[[name]])) return(list(list(valid = TRUE,
      family = object$defaults[[name]]$distribution, distribution = object$defaults[[name]], rank = 0L)))
    Filter(function(c) isTRUE(c$valid), object$candidates[[name]]$alternatives)
  })
  names(options) = names(object$selected)
  fixed = slopes[lengths(object$beliefs[slopes]) == 0L]
  fixed_families = unique(vapply(object$defaults[fixed], `[[`, character(1), "distribution"))
  families = if (length(fixed_families)) fixed_families else c("normal", "student_t")
  # Prefer a complete combination; otherwise retain the largest compatible
  # subset. Sum of mathematical ranks breaks ties, then registry family order.
  coverage = cost = numeric(length(families))
  for (i in seq_along(families)) {
    for (name in slopes) {
      compatible = Filter(function(c) c$family == families[i], options[[name]])
      if (length(compatible)) {
        coverage[i] = coverage[i] + 1L
        cost[i] = cost[i] + compatible[[1]]$rank
      }
    }
  }
  family = if (length(families)) families[order(-coverage, cost)[1]] else NULL
  for (name in names(options)) {
    if (!length(object$beliefs[[name]])) {
      object$selected[name] = list(object$defaults[[name]])
      next
    }
    result = object$candidates[[name]]
    eligible = options[[name]]
    if (name %in% slopes) eligible = Filter(function(c) identical(c$family, family), eligible)
    object$selected[name] = list(NULL)
    if (length(eligible)) {
      chosen = eligible[[which.min(vapply(eligible, `[[`, numeric(1), "rank"))]]
      # Remove diagnostics left over from a previous automatic selection.
      fields = c("family", "distribution", "achieved", "residuals", "score", "rank", "reference_used", "numerical", "convergence")
      result[fields] = lapply(fields, function(field) chosen[[field]])
      result$valid = TRUE
      result$status = "selected"
      result$message = NULL
      object$selected[name] = list(chosen$distribution)
    } else if (length(options[[name]])) {
      result$status = "backend_incompatible"
      result$message = "Valid mathematical candidates exist, but rstanarm requires one slope prior family. Backend defaults were preserved."
    }
    object$candidates[[name]] = result
  }
  object$backend$compatible = all(vapply(object$selected, Negate(is.null), logical(1))) &&
    length(unique(vapply(Filter(Negate(is.null), object$selected[slopes]), `[[`, character(1), "distribution"))) <= 1L
  object
}
