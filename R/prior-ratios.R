prior_ratio_scale = function(scale) !is.null(scale) && scale %in% c("odds_ratio", "rate_ratio")

prior_ratio_validate = function(record, object, parameter) {
  family = object$model$family$family

  # Match the multiplicative interpretation to the model link.
  expected = switch(family, binomial = "odds_ratio", poisson = "rate_ratio", "unsupported")
  if (record$scale != expected) stop("Requested ratio scale does not match this model family: ", record$scale)
  metadata = prior_numeric_slope_metadata(object$model)
  if (parameter == "(Intercept)" || !parameter %in% names(metadata)) {
    stop(
      "Ratio beliefs require an untransformed numeric slope without interactions; unsupported parameter: ",
      parameter
    )
  }
  predictor = metadata[[parameter]]$predictor
  expressions = Filter(is.call, as.list(attr(object$model$terms, "variables"))[-1])
  if (any(vapply(expressions, prior_expression_uses_predictor, logical(1), predictor = predictor))) {
    stop(
      "Ratio beliefs require an untransformed numeric slope without other transformations of the same predictor: ",
      parameter
    )
  }
  value = if (record$type == "probability") record$threshold else record$target
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) stop("Ratio targets and thresholds must be positive finite numbers.")
  record
}

# One exp(normal) constraint can be resolved analytically without an SD cap.
prior_fit_single_ratio = function(record, reference) {
  location = reference$mean
  spread = reference$sd
  retained = paste0("Backend-default link SD = ", format(spread), " retained unchanged; location elicited.")
  if (record$type == "mode") location = log(record$target) + spread^2
  if (record$type == "mean") location = log(record$target) - spread^2 / 2
  if (record$type == "sd") {
    # sd^2 / exp(2*mu) = w*(w-1), w=exp(s^2). Stable at large
    # or small targets, without inventing a location or capping the spread.
    logA = 2 * (log(record$target) - location)
    variance = if (logA > 350) .5 * logA + log(.5 * exp(-.5 * logA) + sqrt(1 + .25 * exp(-logA))) else log1p(2 * exp(logA) / (1 + sqrt(1 + 4 * exp(logA))))
    spread = sqrt(variance)
    retained = paste0("Backend-default link mean = ", format(location), " retained unchanged; SD elicited.")
  }
  prior = list(distribution = "normal", mean = location, sd = spread)
  achieved = tryCatch(prior_distribution_target(record, prior), error = function(e) e)
  if (inherits(achieved, "error")) {
    return(list(
      valid = FALSE,
      message = paste("Single ratio belief is outside reliable numerical range; add another belief.", conditionMessage(achieved))
    ))
  }
  residual = (achieved - record$target) / record$tolerance
  valid = is.finite(residual) && abs(residual) <= 1
  list(
    valid = valid,
    distribution = prior,
    achieved = achieved,
    residuals = residual,
    score = residual^2,
    reference_used = TRUE,
    numerical = paste(retained, "Analytic lognormal evaluation; no numerical integration."),
    message = if (!valid) "Single ratio constraint could not be satisfied at the backend reference; add another belief." else NULL
  )
}

prior_expression_uses_predictor = function(expression, predictor) predictor %in% all.vars(expression)

# Display positive quantities from their logarithms, avoiding overflow and
# avoiding printing numerical underflow as a literal zero ratio.
prior_ratio_format_log = function(logvalue) {
  if (!is.finite(logvalue)) return("not numerically representable")
  if (logvalue > log(.Machine$double.xmax)) return(paste0("exceeds numeric range (log = ", format(logvalue, digits = 4), ")"))
  value = exp(logvalue)
  if (value == 0) return(paste0("below numeric range (log = ", format(logvalue, digits = 4), ")"))
  if (value >= 1e6 || value < 1e-3) return(formatC(value, format = "e", digits = 3))
  format(signif(value, 5), scientific = FALSE, trim = TRUE, big.mark = ",")
}

prior_ratio_diagnostics = function(object, parameter) {
  prior = object$selected[[parameter]]
  if (is.null(prior)) {
    prior_console_line("Induced ratio summary unavailable: no backend-compatible selected prior.")
    return(invisible(NULL))
  }
  quantiles = prior_distribution_operation(prior, "quantile", c(.025, .5, .975))
  formatted = vapply(quantiles, prior_ratio_format_log, character(1))
  prior_console_line(paste0("Selected induced ratio distribution: median = ", formatted[2],
    "; central 95% interval = [", formatted[1], ", ", formatted[3], "]."))
  if (prior$distribution == "normal") {
    prior_console_line(paste0("Ratio mean = ", prior_ratio_format_log(prior$mean + prior$sd^2 / 2),
      "; ratio mode = ", prior_ratio_format_log(prior$mean - prior$sd^2), "."))
  } else {
    prior_console_line("Ratio mean: no finite mean (infinite expectation for exponentiated Student's t).")
    prior_console_line("Ratio mode: no attained positive global mode; the density diverges toward zero.")
  }
  records = object$beliefs[[parameter]]
  keys = vapply(records, prior_outcome_constraint_key, character(1))
  unique = records[!duplicated(keys)]
  candidate = object$candidates[[parameter]]
  if (!isTRUE(candidate$reference_used)) {
    prior_console_line("User beliefs jointly determine the link location and spread; no backend-default parameter was retained as a reference.")
    return(invisible(NULL))
  }
  if (length(unique) == 1L && prior_ratio_scale(unique[[1]]$scale) && unique[[1]]$type != "probability") {
    type = unique[[1]]$type
    if (type == "sd") {
      prior_console_line(paste0(
        "User-supplied constraint: ratio SD. Link spread derived from that constraint; backend-default link location retained = ",
        format(object$defaults[[parameter]]$mean, digits = 5),
        "."
      ))
    } else {
      prior_console_line(paste0(
        "User-supplied constraint: ratio ",
        type,
        ". Link location derived from that constraint; backend-default link SD retained = ",
        format(object$defaults[[parameter]]$sd, digits = 5),
        "."
      ))
    }
  } else {
    location_only = all(vapply(records, prior_ratio_location_constraint, logical(1)))
    if (location_only) {
      prior_console_line(paste0("User beliefs constrain link location; backend-default link SD retained = ", format(object$defaults[[parameter]]$sd, digits = 5), "."))
    } else {
      prior_console_line("User quantile beliefs constrain location and spread jointly; neither link parameter was fixed. Both were resolved using the backend-default reference, not additional user-supplied parameters.")
    }
  }
  invisible(NULL)
}

prior_ratio_location_constraint = function(record) {
  if (record$type == "probability") return(record$target == .5)
  !prior_outcome_scale(record$scale) && record$type %in% c("mean", "mode")
}
