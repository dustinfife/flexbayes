# Provenance is derived from recorded reconciliation rows; it is never fitted
# or stored separately. A shared statement is counted once per participating
# parameter, not once per model-matrix entry or as independent information.
prior_report_changed = function(before, after, compromised = FALSE) {
  if (compromised && !identical(before, after)) return(TRUE)
  if (identical(before$distribution, "normal") && identical(after$distribution, "normal") &&
      is.finite(before$sd) && before$sd > 0) {
    # Compare locations in prior-SD units so predictor units cannot change the
    # source label. This is display precision, never a fitting tolerance.
    difference = c((after$mean - before$mean) / before$sd, after$sd / before$sd - 1)
    return(!isTRUE(all.equal(c(0, 0), unname(difference))))
  }
  !isTRUE(all.equal(before, after))
}

prior_indirect_report = function(object, verified = NULL) {
  parameters = names(object$selected)
  result = list(count = stats::setNames(integer(length(parameters)), parameters),
    changed = stats::setNames(rep(FALSE, length(parameters)), parameters),
    specifications = stats::setNames(rep(list(integer()), length(parameters)), parameters),
    single_parameter = logical(), types = character(), total = 0L, verified = FALSE)
  joint = object$joint
  if (is.null(joint)) return(result)
  X = joint$design
  if (!is.matrix(X) || nrow(X) != length(object$outcomes) || is.null(colnames(X)) || any(!is.finite(X))) return(result)
  details = joint[["coefficient"]]
  sizes = lengths(lapply(object$outcomes, `[[`, "beliefs"))
  result$total = sum(sizes)
  result$types = vapply(object$outcomes, function(s) if (is.null(s$difference)) "outcome" else "difference", character(1))
  result$single_parameter = rowSums(X != 0) == 1L
  for (name in intersect(parameters, colnames(X))) {
    indices = which(X[, name] != 0)
    result$specifications[[name]] = indices
    result$count[name] = sum(sizes[indices])
    ids = vapply(object$beliefs[[name]], `[[`, numeric(1), "id")
    compromised = is.data.frame(details) && any(details$satisfied[match(ids, details$id)] %in% FALSE)
    if (!is.null(joint$baseline[[name]])) result$changed[name] = prior_report_changed(
      joint$baseline[[name]], object$selected[[name]], compromised)
  }
  if (!is.null(verified)) { result$verified = isTRUE(verified); return(result) }
  if (isTRUE(joint$retained) || !joint$status %in% c("compatible", "reconciled")) return(result)
  # Printing reports recorded provenance, without reevaluating formulas or
  # beliefs. Executable validation and ordinary-summary readiness are unchanged.
  result$verified = identical(object$beliefs, joint$beliefs) &&
    identical(object$outcomes, joint$specifications) &&
    identical(object$selected, joint$adjusted) &&
    identical(X, joint$problem$X)
  result
}

prior_report_source = function(parameter, object, indirect) {
  retained = object$fallbacks[[parameter]]
  if (!is.null(retained)) {
    if (retained$source != "backend_default") return("Previous elicited (fallback)")
    return(if (object$backend$package == "brms") "flexbayes default (fallback)" else "Backend default (fallback)")
  }
  candidate = object$candidates[[parameter]]
  if (identical(candidate$status, "contradictory")) return("Contradictory")
  if (identical(candidate$status, "backend_incompatible")) return("Incompatible")
  if (is.null(object$selected[[parameter]])) return("Unresolved")
  if (isTRUE(indirect$changed[parameter])) {
    if (!indirect$verified) return("Unverified")
    indices = indirect$specifications[[parameter]]
    if (!length(object$beliefs[[parameter]]) && length(indices) == 1L && indirect$single_parameter[indices]) {
      return(paste("Elicited via", indirect$types[indices]))
    }
    return("Reconciled")
  }
  if (!length(object$beliefs[[parameter]])) return(if (object$backend$package == "brms") "flexbayes default" else "Backend default")
  if (isTRUE(candidate$reference_used)) return("Elicited*")
  "Elicited"
}

prior_report_beliefs = function(direct, indirect, verified) {
  if (!indirect) return(as.character(direct))
  label = paste(indirect, if (verified) "indirect" else "pending")
  if (!direct) return(label)
  paste(direct, "+", indirect)
}

prior_report_table = function(object, indirect = prior_indirect_report(object)) {
  parameters = names(object$selected)
  counts = lengths(object$beliefs)
  beliefs = if (any(indirect$count)) vapply(parameters, function(name)
    prior_report_beliefs(counts[name], indirect$count[name], indirect$verified), character(1)) else counts
  data.frame(Parameter = parameters, Prior = vapply(object$selected, prior_format, character(1)),
    Source = vapply(parameters, prior_report_source, character(1), object = object, indirect = indirect),
    Beliefs = unname(beliefs), row.names = NULL)
}

prior_report_indirect_note = function(indirect) {
  if (!any(indirect$count)) return(invisible(NULL))
  prefix = if (indirect$verified) "Indirect" else "Pending indirect"
  prior_console_line(paste0(prefix, ": ", indirect$total,
    " statement(s) total. Mixed counts are direct + indirect. Counts show participation, not change or independent information; do not sum shared rows."))
  invisible(NULL)
}

prior_report_unchanged_participants = function(indirect) {
  parameters = names(indirect$count)[indirect$count > 0 & !indirect$changed]
  if (length(parameters)) prior_console_line(paste("Unchanged participants:", paste(parameters, collapse = ", ")))
  invisible(NULL)
}
