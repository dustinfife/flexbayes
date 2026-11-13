prior_fit_candidates = function(records, reference) {
  families = c("normal", "student_t")
  boundary = any(vapply(records, function(r) r$type == "probability" && r$target %in% c(0, 1), logical(1)))
  alternatives = stats::setNames(lapply(families, function(family) {
    candidate = if (boundary) list(valid = FALSE, message = "Boundary probability is consistent, but no supported continuous family attains it at a finite threshold.") else prior_fit_location_scale(records, reference, family)
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

prior_select_compatible = function(object, parameter) {
  result = object$candidates[[parameter]]
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  others = setdiff(slopes, parameter)
  fixed = object$selected[others]
  families = unique(vapply(Filter(Negate(is.null), fixed), `[[`, character(1), "distribution"))
  alternatives = Filter(function(c) isTRUE(c$valid), result$alternatives)
  if (parameter %in% slopes) alternatives = Filter(function(c) !length(families) || (length(families) == 1L && c$family == families), alternatives)
  if (length(alternatives)) {
    chosen = alternatives[[1]]
    result[names(chosen)] = chosen
    result$status = "selected"
    result$message = NULL
    object$selected[parameter] = list(chosen$distribution)
  } else if (isTRUE(result$valid)) {
    result$status = "backend_incompatible"
    result$message = "Valid mathematical candidates exist, but rstanarm requires one slope prior family. Other parameters' existing priors were preserved."
  }
  object$candidates[[parameter]] = result
  object$backend$compatible = all(vapply(object$selected, Negate(is.null), logical(1))) && length(unique(vapply(Filter(Negate(is.null), object$selected[slopes]), `[[`, character(1), "distribution"))) <= 1L
  object
}
