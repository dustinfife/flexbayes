#' Simulate Gaussian prior predictive responses
#'
#' Uses selected priors without posterior fitting. Backend defaults may depend
#' on observed response summaries; their resolved values are kept fixed even
#' with newdata. Parameter draws retain the centered intercept, while
#' coefficient_draws contains the equivalent uncentered intercept.
#' @param object A flexbayes_priors object.
#' @param n Number of simulations.
#' @param newdata Optional predictor data frame; no response is required.
#' @return A flexbayes_prior_predictive object with draws-by-observations
#'   matrices expected and simulated, parameter_draws, coefficient_draws,
#'   the prediction design matrix, and metadata.
#' @export
check_prior = function(object, n = 1000, newdata = NULL) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n) || n > .Machine$integer.max) stop("n must be a positive integer.")
  if (object$model$family$family != "gaussian" || object$model$family$link != "identity") stop("Only Gaussian identity prior predictive checks are supported.")
  prior_predictive_validate(object)
  X = prior_prediction_matrix(object, newdata)
  draws = prior_sample_parameters(object, n)
  coefficients = draws[, colnames(X), drop = FALSE]
  if ("(Intercept)" %in% colnames(X)) {
    slopes = setdiff(colnames(X), "(Intercept)")
    coefficients[, "(Intercept)"] = draws[, "(Intercept)"] -
      as.vector(draws[, slopes, drop = FALSE] %*% object$model$reference[slopes])
  }
  expected = coefficients %*% t(X)
  simulated = prior_simulate_gaussian(expected, draws[, "sigma"])
  if (any(!is.finite(expected)) || any(!is.finite(simulated))) stop("Non-finite prior predictions; selected priors or predictor values exceed numerical limits.")
  result = list(parameter_draws = draws, coefficient_draws = coefficients,
    expected = expected, simulated = simulated, X = X,
    metadata = list(n = n, observations = nrow(X), family = object$model$family,
      formula = object$model$formula, reference = object$model$reference,
      priors = object$selected, backend = object$backend,
      newdata = !is.null(newdata),
      default_scaling = "Backend defaults resolved from original mean(y), sd(y), and design; no posterior conditioning."))
  class(result) = "flexbayes_prior_predictive"
  result
}

prior_predictive_validate = function(object) {
  prior_validate_standardized_beliefs(object)
  required = c(colnames(object$model$X), "sigma")
  if (!all(required %in% names(object$selected))) stop("Required prior is missing.")
  for (name in required) {
    status = object$candidates[[name]]$status
    if (!is.null(status) && status %in% c("contradictory", "backend_incompatible", "no_supported_fit")) stop("Cannot sample prior for ", name, ": ", status, ".")
    prior = object$selected[[name]]
    if (is.null(prior)) stop("No selected prior for ", name, ".")
    tryCatch(prior_validate_distribution(prior), error = function(error) {
      stop("Cannot sample prior for ", name, ": ", conditionMessage(error), call. = FALSE)
    })
    entry = prior_distribution_registry()[[prior$distribution]]
    if (!is.function(entry$random)) stop("Prior cannot be sampled: ", name)
    if (name == "sigma" && entry$support[1] < 0) stop("Residual-scale prior must have positive support.")
  }
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  families = unique(vapply(object$selected[slopes], `[[`, character(1), "distribution"))
  if (!isTRUE(object$backend$compatible) || length(families) > 1L) stop("Prior specification is incompatible with the backend.")
  if ("(Intercept)" %in% required && !object$selected$`(Intercept)`$distribution %in% c("normal", "student_t")) stop("Unsupported backend intercept prior.")
  if (object$backend$package != "rstanarm" || any(!families %in% c("normal", "student_t"))) stop("Unsupported backend prior specification.")
  invisible(object)
}

prior_prediction_matrix = function(object, newdata) {
  if (is.null(newdata)) return(object$model$X)
  if (!is.data.frame(newdata) || !nrow(newdata)) stop("newdata must be a nonempty data frame.")
  terms = stats::delete.response(object$model$terms)
  variables = all.vars(terms)
  if (!all(variables %in% names(newdata))) stop("newdata is missing predictors: ", paste(setdiff(variables, names(newdata)), collapse = ", "))
  factors = vapply(object$model$frame, is.factor, logical(1))
  levels = lapply(object$model$frame[factors], levels)
  frame = stats::model.frame(terms, newdata, na.action = stats::na.fail, xlev = levels)
  X = stats::model.matrix(terms, frame, contrasts.arg = object$model$contrasts)
  if (!identical(colnames(X), colnames(object$model$X)) || nrow(X) != nrow(newdata) || any(!is.finite(X))) stop("newdata does not produce a finite compatible design matrix.")
  X
}

prior_sample_parameters = function(object, n) {
  required = c(colnames(object$model$X), "sigma")
  draws = matrix(NA_real_, nrow = n, ncol = length(required), dimnames = list(NULL, required))
  for (name in required) {
    values = tryCatch(prior_distribution_operation(object$selected[[name]], "random", n),
      error = function(error) stop("Cannot sample prior for ", name, ": ", conditionMessage(error), call. = FALSE))
    if (length(values) != n || any(!is.finite(values))) stop("Non-finite or invalid sampled values for ", name)
    draws[, name] = values
  }
  if (any(draws[, "sigma"] <= 0)) stop("Residual-scale draws must be positive.")
  draws
}

prior_simulate_gaussian = function(expected, sigma) {
  noise = matrix(stats::rnorm(length(expected)), nrow = nrow(expected))
  expected + noise * sigma
}

#' Plot prior predictive intervals
#' @param x A flexbayes_prior_predictive object.
#' @param type Show both expected-response uncertainty and individual-observation
#'   variability (default), or only expected or simulated responses.
#' @param ... Reserved.
#' @return A ggplot object showing medians and central 90 percent intervals by
#'   prediction row. Row order does not imply a continuous predictor.
#' @export
plot.flexbayes_prior_predictive = function(x, type = c("both", "expected", "simulated"), ...) {
  type = match.arg(type)
  types = if (type == "both") c("expected", "simulated") else type
  rows = list()
  for (kind in types) {
    values = x[[kind]]
    intervals = apply(values, 2, stats::quantile, probs = c(.05, .5, .95), names = FALSE)
    rows[[kind]] = data.frame(row = seq_len(ncol(values)), lower = intervals[1, ],
      median = intervals[2, ], upper = intervals[3, ],
      response = if (kind == "expected") "Expected response" else "Individual observation")
  }
  ggplot2::ggplot(do.call(rbind, rows), ggplot2::aes(x = row, y = median)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), alpha = .5) +
    ggplot2::geom_point() + ggplot2::facet_wrap(~response) +
    ggplot2::labs(x = "Prediction row", y = "Response", subtitle = "Prior predictive medians and central 90% intervals")
}
