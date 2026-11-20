prior_fit_candidates = function(records, reference) {
  if (reference$distribution == "exponential") return(prior_fit_scale_candidates(records, reference))
  # Fit each supported family against the complete collection of beliefs.
  families = c("normal", "student_t")
  boundary = any(vapply(records, function(r) r$type == "probability" && r$target %in% c(0, 1), logical(1)))
  alternatives = stats::setNames(lapply(families, function(family) {
    outcome = any(vapply(records, function(r) prior_outcome_scale(r$scale), logical(1)))
    candidate = if (outcome) prior_fit_outcome(records, reference, family) else if (boundary) list(
      valid = FALSE,
      message = "Boundary probability is consistent, but no supported continuous family attains it at a finite threshold."
    ) else prior_fit_location_scale(records, reference, family)
    candidate$family = family
    candidate
  }), families)

  valid = which(vapply(alternatives, function(c) isTRUE(c$valid), logical(1)))

  # Once tolerances pass, prefer the simpler normal model. Fixed-df t has
  # equal fitted dimension but adds a tail-shape assumption; numerical noise
  # in near-zero fit scores must not decide the family.
  for (i in seq_along(valid)) alternatives[[valid[i]]]$rank = i
  best = if (length(valid)) alternatives[[valid[1]]] else list(
    valid = FALSE,
    message = paste("No supported distribution fits all beliefs within tolerance.", paste(unique(vapply(alternatives, `[[`, character(1), "message")), collapse = " "))
  )
  best$status = if (length(valid)) "selected" else "no_supported_fit"
  best$alternatives = alternatives
  best$mathematical = if (length(valid)) alternatives[[valid[1]]]$distribution else NULL
  best
}

# Recompute automatic selections without refitting or changing candidate ranks.
# Defaults are fixed; elicited candidates are automatic choices, not user locks.
prior_choose_candidates = function(object) {
  options = lapply(names(object$selected), function(name) {
    if (!length(object$beliefs[[name]])) return(list(list(valid = TRUE,
      family = object$defaults[[name]]$distribution, distribution = object$defaults[[name]], rank = 0L)))
    Filter(function(c) isTRUE(c$valid), object$candidates[[name]]$alternatives)
  })
  names(options) = names(object$selected)
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  families = c("normal", "student_t")
  coverage = cost = numeric(length(families))
  feasible = rep(TRUE, length(families))

  # A family must accommodate every actual prior, including retained fallbacks.
  for (i in seq_along(families)) {
    for (name in slopes) {
      eligible = Filter(function(c) c$family == families[i], options[[name]])
      if (length(eligible)) {
        coverage[i] = coverage[i] + 1L
        cost[i] = cost[i] + eligible[[1]]$rank
      } else if (object$last_valid[[name]]$distribution != families[i]) feasible[i] = FALSE
    }
  }
  available = which(feasible)
  if (!length(available)) stop("No coherent retained slope-prior collection is available; identify priors again.")
  family = families[available[order(-coverage[available], cost[available])[1]]]
  object$fallbacks = list()

  for (name in names(options)) {
    result = object$candidates[[name]]
    eligible = options[[name]]
    if (name %in% slopes) eligible = Filter(function(c) c$family == family, eligible)
    if (length(eligible)) {
      chosen = eligible[[which.min(vapply(eligible, `[[`, numeric(1), "rank"))]]
      object$selected[name] = list(chosen$distribution)
      if (!length(object$beliefs[[name]])) {
        object$candidates[[name]] = list()
        next
      }
      fields = c("family", "distribution", "achieved", "residuals", "score", "rank", "reference_used", "numerical", "convergence")
      result[fields] = lapply(fields, function(field) chosen[[field]])
      result$valid = TRUE
      result$status = "selected"
      result$message = NULL
    } else {
      if (length(options[[name]])) {
        result$status = "backend_incompatible"
        result$message = "Valid mathematical candidates exist, but rstanarm requires one slope prior family. A previous compatible prior has been retained."
      }
      if (is.null(result$status)) result$status = "no_supported_fit"
      object$selected[name] = list(object$last_valid[[name]])
      object$fallbacks[[name]] = prior_fallback_record(object, name, result$status, result$message)
    }
    object$candidates[[name]] = result
  }
  object$backend$compatible = prior_available_compatible(object)
  object
}
