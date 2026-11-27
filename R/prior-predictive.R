#' Simulate prior predictive responses
#'
#' Uses selected priors without posterior fitting. Backend defaults may depend
#' on observed response summaries; their resolved values are kept fixed even
#' with newdata. Parameter draws retain the centered intercept, while
#' coefficient_draws contains the equivalent uncentered intercept. Expected
#' responses use the identity, logistic, or exponential inverse link. Poisson
#' simulation rejects nonfinite rates or rates above 1e8.
#' Reconciled outcome priors use the same simulation machinery. Use
#' response='expected' to examine beliefs without individual residual noise;
#' simulated outcomes additionally use the selected residual-scale prior.
#' Expected-outcome updates never change sigma; update_priors(parameter = 'sigma',
#' ...) can elicit it separately.
#' @param object A flexbayes_priors object.
#' @param allow_fallback Use retained priors despite unsuccessful elicitation?
#'   Defaults to FALSE. TRUE explicitly approves their use and issues a warning;
#'   current beliefs remain recorded and are not treated as satisfied. Model-data
#'   validation and backend restrictions still apply.
#' @param n Number of simulations.
#' @param newdata Optional predictor data frame; no response is required.
#'   Simulations use these rows, or the original estimation cases when omitted.
#'   Model plots use a separate balanced grid unless newdata is supplied.
#' @return A flexbayes_prior_predictive object with draws-by-observations
#'   matrices expected and simulated, parameter_draws, coefficient_draws,
#'   the linear_predictor and prediction design matrix, and metadata retaining
#'   prediction data, frozen model terms, factor levels and contrasts for plotting,
#'   and a compact deterministic model-grid recipe. Creating the recipe does not
#'   change simulation cases, parameter draws, or the random-number state.
#' @export
check_prior = function(object, n = 1000, newdata = NULL, allow_fallback = FALSE) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n) || n > .Machine$integer.max) stop("n must be a positive integer.")
  prior_validate_family(object$model$family)
  prior_predictive_validate(object, allow_fallback)

  # Sample the selected priors once, keeping joint draw identities intact.
  X = prior_prediction_matrix(object, newdata)
  draws = prior_sample_parameters(object, n)
  coefficients = draws[, colnames(X), drop = FALSE]

  # Convert the centered intercept only for prediction on the original data.
  if ("(Intercept)" %in% colnames(X)) {
    slopes = setdiff(colnames(X), "(Intercept)")
    coefficients[, "(Intercept)"] = draws[, "(Intercept)"] -
      as.vector(draws[, slopes, drop = FALSE] %*% object$model$reference[slopes])
  }

  # Expected responses and individual observations share the same parameters.
  linear = coefficients %*% t(X)
  response = prior_simulate_response(linear, draws, object$model$family)
  expected = response$expected
  simulated = response$simulated
  if (any(!is.finite(expected)) || any(!is.finite(simulated))) stop("Non-finite prior predictions; selected priors or predictor values exceed numerical limits.")

  # Save the original model information for later plotting and validation.
  result = list(
    parameter_draws = draws,
    coefficient_draws = coefficients,
    expected = expected,
    simulated = simulated,
    linear_predictor = linear,
    X = X,
    metadata = list(n = n, observations = nrow(X), family = object$model$family,
      formula = object$model$formula, reference = object$model$reference,
      priors = object$selected, backend = object$backend,
      newdata = !is.null(newdata),
      prediction_data = if (is.null(newdata)) object$model$data else as.data.frame(newdata),
      observed_response = stats::model.response(object$model$frame),
      terms = object$model$terms, contrasts = object$model$contrasts,
      factor_levels = stats::.getXlevels(object$model$terms, object$model$frame),
      default_scaling = prior_default_description(object))
  )
  result$metadata$model_grid = prior_model_grid_specification(object, newdata)
  class(result) = "flexbayes_prior_predictive"
  if (length(prior_fallback_issues(object))) result$metadata$fallback = list(
    approved = allow_fallback, parameters = object$fallbacks, outcome_status = object$joint$status)
  result
}

prior_predictive_validate = function(object, allow_fallback = FALSE) {
  issues = prior_fallback_permission(object, allow_fallback)
  if (length(issues)) prior_validate_retained(object)
  if (!isTRUE(object$joint$retained)) prior_joint_validate(object)
  prior_r_squared_validate(object, allow_fallback)
  prior_validate_standardized_beliefs(object)
  required = prior_required_parameters(object)
  if (!all(required %in% names(object$selected))) stop("Required prior is missing.")
  for (name in required) {
    status = object$candidates[[name]]$status
    if (!is.null(status) && status %in% c("contradictory", "backend_incompatible", "no_supported_fit") &&
        !(allow_fallback && !is.null(object$fallbacks[[name]]))) stop("Cannot sample prior for ", name, ": ", status, ".")
    prior = object$selected[[name]]
    if (is.null(prior)) stop("No selected prior for ", name, ".")
    tryCatch(prior_validate_distribution(prior), error = function(error) {
      stop("Cannot sample prior for ", name, ": ", conditionMessage(error), call. = FALSE)
    })
    entry = prior_distribution_registry()[[prior$distribution]]
    if (!is.function(entry$random)) stop("Prior cannot be sampled: ", name)
    if (name == "sigma" && prior_has_sigma(object) && entry$support[1] < 0) stop("Residual-scale prior must have positive support.")
  }
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  families = unique(vapply(object$selected[slopes], `[[`, character(1), "distribution"))
  if (!isTRUE(object$backend$compatible) || (object$backend$package == "rstanarm" && length(families) > 1L)) stop("Prior specification is incompatible with the backend.")
  if ("(Intercept)" %in% required && !object$selected$`(Intercept)`$distribution %in% c("normal", "student_t")) stop("Unsupported backend intercept prior.")
  if (!object$backend$package %in% c("rstanarm", "brms") || any(!families %in% c("normal", "student_t"))) stop("Unsupported backend prior specification.")
  if (prior_has_sigma(object) && !object$selected$sigma$distribution %in% prior_scale_families(object$backend$package)) stop("Unsupported sigma prior for selected backend.")
  if (object$backend$package == "brms") prior_brms_validate_data(object)
  if (length(issues)) prior_warn_fallback(object)
  invisible(object)
}

prior_prediction_matrix = function(object, newdata) {
  if (is.null(newdata)) return(object$model$X)
  if (!is.data.frame(newdata) || !nrow(newdata)) stop("newdata must be a nonempty data frame.")
  terms = stats::delete.response(object$model$terms)
  variables = all.vars(terms)
  if (!all(variables %in% names(newdata))) stop("newdata is missing predictors: ", paste(setdiff(variables, names(newdata)), collapse = ", "))

  # Retain the estimation sample's factor levels and contrasts for new cases.
  levels = stats::.getXlevels(object$model$terms, object$model$frame)
  frame = stats::model.frame(terms, newdata, na.action = stats::na.fail, xlev = levels)
  X = stats::model.matrix(terms, frame, contrasts.arg = object$model$contrasts)
  if (!identical(colnames(X), colnames(object$model$X)) || nrow(X) != nrow(newdata) || any(!is.finite(X))) stop("newdata does not produce a finite compatible design matrix.")
  X
}

prior_sample_parameters = function(object, n) {
  required = prior_required_parameters(object)
  draws = matrix(NA_real_, nrow = n, ncol = length(required), dimnames = list(NULL, required))
  for (name in required) {
    values = tryCatch(prior_distribution_operation(object$selected[[name]], "random", n),
      error = function(error) stop("Cannot sample prior for ", name, ": ", conditionMessage(error), call. = FALSE))
    if (length(values) != n || any(!is.finite(values))) stop("Non-finite or invalid sampled values for ", name)
    draws[, name] = values
  }
  if (prior_has_sigma(object) && any(draws[, "sigma"] <= 0)) stop("Residual-scale draws must be positive.")
  draws
}

prior_simulate_gaussian = function(expected, sigma) {
  noise = matrix(stats::rnorm(length(expected)), nrow = nrow(expected))
  expected + noise * sigma
}

#' Plot prior predictive distributions and relationships
#'
#' Plot outcome distributions, expected relationships, or per-case intervals
#' from `check_prior()`. For individual parameter densities, use `plot(priors)`
#' or `plot(priors, parameter = "x")` instead.
#' @param x A flexbayes_prior_predictive object.
#' @param type Outcome distribution (default), model relationships, or interval
#'   summaries: outcome, model, or interval.
#' @param response simulated (default for outcome), expected (default for model),
#'   or both (interval only; its default).
#' @param formula Optional literal flexplot visualization formula for model plots.
#' @param draws Maximum individual joint prior draws shown in draws/both model
#'   displays (default 30); never limits draws used to calculate intervals.
#' @param display Model display: interval (default), draws, or both. Intervals
#'   summarize prior uncertainty in the expected response, not individual outcomes.
#' @param max_display Display budget for outcome values or model-grid rows.
#'   Outcome plots select evenly spaced matrix indices deterministically, retaining
#'   draw/case IDs; simulated model plots reduce draw count to fit all cases.
#'   Stored simulations are never changed. Interval plots retain all cases.
#' @param prediction.alpha Opacity of prior relationships and simulated points.
#' @param bins,breaks,labels Model-display binning controls passed to prepare_flexplot.
#' @param ... Reserved; unsupported arguments are rejected for new plot types.
#' @details Outcome plots pool across draws and the supplied predictor cases,
#'   not an arbitrary population. Continuous values use 30-bin histograms; simulated
#'   binary/count values use discrete proportion bars. All plots use theme_bw().
#'   Captions count joint parameter draws, not pooled outcome values. Model plots use public
#'   flexplot preparation contract 1 and prediction.draw rendering. Summary
#'   displays calculate pointwise medians and 5th/95th percentiles using all stored
#'   joint coefficient draws, after the inverse link on the response scale. Numeric
#'   x uses a median line and translucent ribbon; categorical x uses median points
#'   and vertical intervals. Both adds the selected individual relationships.
#'   Individual relationships use lines for numeric x and unconnected points
#'   for categorical x, retaining flexplot's native positions and draw identities.
#'   Intervals are pointwise prior intervals for expected responses, not
#'   simultaneous bands, posterior intervals, or individual-observation intervals.
#'   Plot-ready layers reuse flexplot's returned prediction mappings, binned data,
#'   and position objects; no binning or dodge logic is reconstructed. This small
#'   dependency on returned layer structure is guarded and covered by rendering
#'   tests. An unsupported layer structure errors and suggests display='draws'.
#'   Original contrasts, levels and transformation terms use stored coefficient
#'   draws; no model is fitted or new parameters/outcomes sampled.
#'
#'   Model plots add subtle dashed lines labelled Observed minimum and Observed
#'   maximum, using the original observed response range even with newdata.
#'   These global references repeat across facets and naturally train the y axis.
#'   Missing/nonfinite observations are ignored; when no original finite response
#'   is available, references are omitted with a warning.
#'   Observed points are hidden with raw.data=FALSE and sample=0 for expected
#'   model plots. flexplot still retains its ordinary axis extent, category-order
#'   conventions and categorical-binary predictor-case counts. Its current public
#'   interface cannot train model y axes exclusively from predictions: hidden
#'   observed responses can still expand the axis range. For predictor-only
#'   newdata, stored expected medians supply the preparer's required outcome
#'   column; no observed response is overwritten. Simulated mode overlays stored
#'   observations at their actual cases with expected marginal curves.
#'
#'   Automatic model grids include every raw model predictor. Continuous nuisance
#'   predictors use five equally spaced values over their observed range; numeric
#'   x axes use up to 51. Integer or low-cardinality predictors use observed unique
#'   values (up to nine for nuisance predictors), and categories retain all levels.
#'   Full combinations receive equal weight, rather than observed-case frequency
#'   weights. Omitted predictors are averaged within each complete coefficient
#'   draw, after applying the inverse link. Medians and central 90% intervals are
#'   then calculated across all draws. Binned groups average their actual grid
#'   values, without predicting at fractional discrete values.
#'
#'   Grids are capped at 8,192 cases. Non-x numeric resolutions are reduced
#'   deterministically, no lower than two, with a warning; x resolution and category
#'   levels are preserved. If the cap still cannot be met, supply explicit newdata.
#'   Explicit newdata is never expanded: omitted predictors are averaged over
#'   supplied rows sharing the displayed coordinates. Inspect the evaluated cases,
#'   resolutions and weighting through attr(plot, 'prior_model_grid'). Unsupported
#'   joint configurations (for example, invalid transformed predictors) require
#'   explicit newdata restricting predictions to valid cases.
#'
#'   Intercept-only models have outcome/interval plots, not model paths. Transformed
#'   responses, transformed visualization formulas, missing preparation rows,
#'   excessive grids and unsupported flexplot configurations fail clearly.
#'   Old saved checks need to be regenerated for model plotting.
#' @return A genuine ggplot object. Plotting does not alter stored simulations.
#' @seealso [plot.flexbayes_priors()], [check_prior()]
#' @examples
#' \dontrun{
#' priors = identify_priors(mpg ~ wt + cyl, data = mtcars)
#' plot(priors)
#' plot(priors, parameter = "wt")
#' check = check_prior(priors, n = 500)
#' plot(check, type = "outcome")
#' plot(check, type = "outcome", response = "expected")
#' plot(check, type = "model") # median and central 90% prior interval
#' plot(check, type = "model", display = "draws")
#' plot(check, type = "model", display = "both")
#' plot(check, type = "model", formula = mpg ~ wt | cyl)
#' plot(check, type = "interval")
#' plot(check, type = "interval", response = "expected")
#' }
#' @export
plot.flexbayes_prior_predictive = function(x, type = c("outcome", "model", "interval"),
    response = NULL, formula = NULL, draws = 30, max_display = 100000,
    prediction.alpha = .15, bins = 3, breaks = NULL, labels = NULL,
    display = c("interval", "draws", "both"), ...) {
  type = match.arg(type)
  if (type == "interval") return(prior_plot_interval(x, if (is.null(response)) "both" else match.arg(response, c("both", "expected", "simulated"))))
  if (is.null(response)) response = if (type == "model") "expected" else "simulated"
  response = match.arg(response, c("expected", "simulated"))
  prior_display_integer(max_display, "max_display")
  if (length(list(...))) stop("Unsupported plotting arguments in ...; use the documented display controls.")
  if (type == "outcome") return(prior_plot_outcome(x, response, max_display))
  prior_plot_model(x, response, formula, draws, max_display, prediction.alpha, bins, breaks, labels, match.arg(display))
}

prior_plot_interval = function(x, response = "both") {
  types = if (response == "both") c("expected", "simulated") else response
  rows = list()
  for (kind in types) {
    values = x[[kind]]
    intervals = apply(values, 2, stats::quantile, probs = c(.05, .5, .95), names = FALSE)
    rows[[kind]] = data.frame(
      row = seq_len(ncol(values)),
      lower = intervals[1, ],
      median = intervals[2, ],
      upper = intervals[3, ],
      response = if (kind == "expected") "Expected outcome" else "Predicted outcome"
    )
  }
  ggplot2::ggplot(do.call(rbind, rows), ggplot2::aes(x = row, y = median)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), alpha = .5) +
    ggplot2::geom_point() + ggplot2::facet_wrap(~response) + ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Prediction row",
      y = switch(x$metadata$family$family, binomial = "Probability / binary response", poisson = "Expected count / count", "Response"),
      title = "Prior predictive intervals"
    )
}
