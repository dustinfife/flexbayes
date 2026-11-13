#' Describe a prior belief
#' @param expression A comparison of value with a numeric constant.
#' @param probability Probability of the comparison.
#' @param mean,mode,sd Numeric moment or mode target; supply one.
#' @param tolerance Absolute acceptance tolerance on the target scale.
#' @return A prior_belief record.
#' @export
belief = function(expression, probability = NULL, mean = NULL, mode = NULL, sd = NULL, tolerance = 1e-6) {
  expr = substitute(expression)
  targets = list(mean = mean, mode = mode, sd = sd)
  supplied = !vapply(targets, is.null, logical(1))
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) stop("Positive finite tolerance required.")
  if (!missing(expression)) {
    if (any(supplied) || !is.call(expr) || length(expr) != 3L || !as.character(expr[[1]]) %in% c("<", ">", "<=", ">=") || !identical(expr[[2]], as.name("value"))) stop("Use value < constant or value > constant.")
    threshold = expr[[3]]
    if (is.call(threshold) && length(threshold) == 2L && identical(threshold[[1]], as.name("-")) && is.numeric(threshold[[2]])) threshold = -threshold[[2]]
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold)) stop("Threshold must be a finite numeric literal.")
    if (!is.numeric(probability) || length(probability) != 1 || !is.finite(probability) || probability < 0 || probability > 1) stop("Probability must be between zero and one.")
    out = list(type = "probability", target = probability, threshold = threshold, operator = as.character(expr[[1]]), expression = expr, label = paste(deparse(expr), "with probability", probability))
  } else {
    if (sum(supplied) != 1L || !is.null(probability)) stop("Supply exactly one mean, mode, or sd.")
    type = names(targets)[supplied]
    target = targets[[type]]
    if (!is.numeric(target) || length(target) != 1 || !is.finite(target) || (type == "sd" && target <= 0)) stop("Invalid belief target.")
    out = list(type = type, target = target, label = paste(type, "=", target))
  }
  out$tolerance = tolerance
  class(out) = "prior_belief"
  out
}

prior_fit_location_scale = function(records, reference, family = "normal") {
  for (type in c("mean", "mode", "sd")) {
    same = Filter(function(r) r$type == type, records)
    if (length(same) > 1L) {
      lower = max(vapply(same, function(r) r$target - r$tolerance, numeric(1)))
      upper = min(vapply(same, function(r) r$target + r$tolerance, numeric(1)))
      if (lower > upper) return(list(valid = FALSE, message = paste("Inconsistent", type, "beliefs: target tolerances do not overlap.")))
    }
  }
  entry = prior_distribution_registry()[[family]]
  unit = if (family == "normal") list(distribution = "normal", mean = 0, sd = 1) else list(distribution = "student_t", location = 0, scale = 1, df = prior_student_t_policy()$df)
  A = matrix(numeric(), ncol = 2)
  b = numeric()
  for (r in records) {
    row = switch(r$type, mean = c(1, 0), mode = c(1, 0), sd = c(0, entry$sd_factor), probability = c(1, prior_distribution_operation(unit, "quantile", if (r$operator %in% c("<", "<=")) r$target else 1 - r$target)))
    A = rbind(A, row)
    b = c(b, if (r$type == "probability") r$threshold else r$target)
  }
  rank = qr(A)$rank
  base = c(reference$mean, reference$sd / entry$sd_factor)
  if (rank == 1L) {
    # Closest standardized location/scale to the backend default.
    direction = A[1, ]
    solution = base + direction * (b[1] - sum(direction * base)) / sum(direction^2)
  } else solution = as.numeric(qr.solve(A, b))
  if (any(!is.finite(solution)) || solution[2] <= 0) return(list(valid = FALSE, reference_used = rank < 2, message = "No location-scale prior with positive scale satisfies these beliefs."))
  distribution = if (family == "normal") list(distribution = family, mean = solution[1], sd = solution[2]) else list(distribution = family, location = solution[1], scale = solution[2], df = prior_student_t_policy()$df)
  achieved = vapply(records, prior_distribution_target, numeric(1), prior = distribution)
  residuals = (achieved - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  valid = all(abs(residuals) <= 1)
  list(valid = valid, reference_used = rank < 2, distribution = distribution, achieved = achieved, residuals = residuals, score = sum(residuals^2), rank = if (valid) 1L else NA_integer_, message = if (!valid) "Beliefs cannot be satisfied by this distribution family within tolerance." else NULL)
}

prior_distribution_target = function(record, prior) {
  if (record$type != "probability") return(prior_distribution_operation(prior, record$type))
  p = prior_distribution_operation(prior, "cdf", record$threshold)
  if (record$operator %in% c("<", "<=")) p else 1 - p
}

prior_fit_normal = function(records, reference) prior_fit_location_scale(records, reference, "normal")

#' Revise prior beliefs
#' @param object Prior specification.
#' @param ... Belief objects.
#' @param parameter Exact parameter name.
#' @param remove Displayed belief numbers to remove before adding beliefs.
#' @return A revised prior specification.
#' @export
update.flexbayes_priors = function(object, ..., parameter, remove = integer()) {
  if (length(parameter) != 1 || !parameter %in% names(object$beliefs)) stop("Unknown parameter.")
  if (parameter == "sigma") stop("Residual-scale elicitation is not supported yet; backend default retained.")
  additions = list(...)
  if (any(!vapply(additions, inherits, logical(1), "prior_belief"))) stop("All additions must be belief objects.")
  records = object$beliefs[[parameter]]
  if (!is.numeric(remove) || anyNA(remove) || any(!is.finite(remove)) || any(remove != floor(remove)) || any(remove < 1 | remove > length(records)) || anyDuplicated(remove)) stop("Invalid removal numbers.")
  if (length(remove)) records = records[-remove]
  for (r in additions) {
    r$id = object$next_id
    object$next_id = object$next_id + 1L
    records = c(records, list(r))
  }
  object$beliefs[[parameter]] = records
  if (!length(records)) {
    object$selected[parameter] = list(object$defaults[[parameter]])
    object$candidates[parameter] = list(list())
    return(prior_select_compatible(object))
  }
  consistency = prior_check_probability_constraints(records)
  if (consistency$contradictory) {
    object$candidates[[parameter]] = list(valid = FALSE, status = "contradictory",
      conflicts = consistency$conflicts, message = consistency$message)
    object$selected[parameter] = list(NULL)
    return(prior_select_compatible(object))
  }
  result = prior_fit_candidates(records, object$defaults[[parameter]])
  object$candidates[[parameter]] = result
  object$selected[parameter] = list(NULL)
  prior_select_compatible(object)
}
