#' Specify properties of a prior distribution
#'
#' Describe what you believe using familiar quantities: a mean, a most likely
#' value (mode), a standard deviation, a median, or a probability. flexbayes
#' uses these statements as constraints when selecting a prior distribution.
#' One belief may leave several distributions possible; combine beliefs to
#' describe both where the distribution belongs and how widely it spreads.
#'
#' @param expression A comparison with `value` on the left, such as
#'   `value < 100`, paired with `probability`. Supported operators are `<`,
#'   `<=`, `>`, and `>=`. Omit this argument when supplying `mean`, `mode`,
#'   or `sd`. The threshold can be a number, variable, or numeric expression.
#' @param probability The probability assigned to the comparison in
#'   `expression`: one finite number between 0 and 1, inclusive.
#' @param mean The desired mean of the distribution, as one finite number.
#' @param mode The desired most likely value of the distribution, as one
#'   finite number. This is not necessarily its mean.
#' @param sd The desired standard deviation of the distribution, as one
#'   positive finite number. This describes spread, not a location.
#' @param tolerance Positive finite absolute numerical acceptance tolerance,
#'   measured in the target's units. The default is `1e-6`. This controls how closely
#'   a candidate must match the stated target; it is not an estimate of the
#'   researcher's uncertainty. See "Matching an answer" below.
#' @param scale The scale on which the belief is expressed. Defaults to
#'   `"original"`. Other choices are `"standardized"`, `"probability"`,
#'   `"response"`, `"odds_ratio"`, and `"rate_ratio"`, with the model and
#'   parameter restrictions described below.
#'
#' @details
#' ## What does a belief say?
#'
#' `belief(mean = 100)` places the mean at 100. `belief(mode = 100)` places
#' the peak there instead. `belief(sd = 20)` describes the width, without
#' saying where the distribution should peak.
#'
#' For probability statements, `value` stands for the uncertain quantity.
#' `belief(value < 150, probability = .95)` assigns 95% probability below
#' 150. It does not say that 95% of people have outcomes below 150 when the
#' belief concerns a regression coefficient.
#'
#' `belief(value < 100, probability = .50)` specifies 100 as a median for
#' a continuous distribution. There is no separate `median` argument.
#' Strict and inclusive comparisons are retained as entered; they need not
#' mean the same thing for a distribution with probability mass at the threshold.
#'
#' Each call supplies either one of `mean`, `mode`, or `sd`, or one comparison
#' together with `probability`. Only the threshold expression is evaluated,
#' in the calling environment; it must give one finite real numeric value.
#' The comparison itself is captured, and its original expression is retained.
#'
#' ## Combining beliefs
#'
#' `belief()` records a statement; it does not fit a distribution by itself.
#' Add statements to an identified prior object with [update_priors()], or
#' use [`update()`][update.flexbayes_priors] to revise existing beliefs. Adding beliefs
#' keeps earlier ones. To remove a belief, use its current displayed number
#' within the parameter's belief list.
#'
#' A mode alone does not determine a distribution's width. When supported,
#' flexbayes uses the backend's default prior as a disclosed reference to
#' resolve missing information. That is different from knowledge supplied by
#' the researcher. Use `summary(priors, diagnostics = TRUE)` to inspect the
#' reference choice and the achieved constraints.
#'
#' More beliefs do not guarantee a fit. Contradictory statements are retained
#' and reported, as are otherwise consistent statements that the supported
#' families or backend cannot represent. A valid belief record does not by
#' itself establish that a usable prior exists.
#'
#' ## Matching an answer
#'
#' A standard deviation or probability belief describes uncertainty about the
#' quantity. `tolerance` instead describes how closely the fitted distribution
#' must match that answer. It is separate from optimizer convergence settings.
#' The absolute default `1e-6` is a numerical matching rule, not inferred
#' precision in an expert's judgment. It is not silently widened.
#'
#' If an expected-price estimate means "about 300,000, to the nearest 10,000,"
#' `belief(mode = 300000, tolerance = 5000)` records that precision explicitly.
#' A probability reported to the nearest percentage point can use
#' `tolerance = .005`. Establish that precision during elicitation, rather
#' than widening tolerances afterwards to make conflicting answers fit.
#' Location and SD tolerances use the entered units; probability tolerances
#' use probability units. Transform location and SD tolerances along with
#' their targets when changing units.
#'
#' ## Available scales
#'
#' \describe{
#'   \item{`"original"`}{The model's coefficient units: response units per
#'   predictor unit for Gaussian slopes, log-odds for logistic coefficients,
#'   and log-rates for Poisson coefficients. Gaussian `sigma` beliefs describe
#'   the residual standard deviation in response units. An intercept refers to the
#'   recorded model-matrix centering values when the model has an intercept.
#'   For example, `belief(mode = 3)` can describe a Gaussian slope.
#'   With `update_outcome()`, this same scale refers to expected-response
#'   units at the specified predictor values instead.}
#'   \item{`"standardized"`}{A fully standardized Gaussian identity-model
#'   slope: `beta_star = beta * sd(x) / sd(y)`. For example,
#'   `belief(mode = .3, scale = "standardized")`. SDs use R's `sd()`
#'   convention and the observed estimation sample. Only plain numeric
#'   slopes without interactions are supported, not intercepts, factor
#'   coefficients, or transformed terms. Predictor and outcome SDs must be
#'   positive and finite. The beliefs are converted to coefficient units;
#'   model data are not standardized.}
#'   \item{`"probability"`}{The success probability at the centered-intercept
#'   reference of a binomial logit model. Only `"(Intercept)"` is supported.
#'   For example, `belief(value < .8, probability = .95,
#'   scale = "probability")` puts 95% prior probability below a success
#'   probability of .8. Mean and mode targets must be strictly between 0
#'   and 1, SD targets strictly between 0 and .5, and comparison thresholds
#'   between 0 and 1, inclusive.}
#'   \item{`"response"`}{The expected count at the centered-intercept
#'   reference of a Poisson log model. Only `"(Intercept)"` is supported.
#'   For example, `belief(mean = 3, scale = "response")` describes the
#'   mean of the distribution of expected counts, not a link-scale mean.
#'   Moment and mode targets must be positive; comparison thresholds may
#'   be zero. This does not describe individual count variability.}
#'   \item{`"odds_ratio"`}{The multiplicative change in odds for a one-unit
#'   increase in a plain numeric slope of a binomial logit model, holding
#'   other predictors fixed. For example, `belief(value < 2,
#'   probability = .95, scale = "odds_ratio")`. Intercepts, factors,
#'   interactions, transformed terms, and predictors also used in another
#'   transformation are unsupported. Targets and thresholds must be positive.}
#'   \item{`"rate_ratio"`}{The multiplicative change in expected count for a
#'   one-unit increase in a plain numeric slope of a Poisson log model.
#'   For example, `belief(mode = 1.5, scale = "rate_ratio")`. The same
#'   parameter restrictions and positive-target requirements as for
#'   `"odds_ratio"` apply.}
#' }
#'
#' Beliefs on different supported scales can be combined for the same parameter.
#' On transformed scales, means, modes, and SDs describe the induced
#' distribution, not simply a transformation of a link-scale moment. Requested
#' moments or modes must exist for the candidate family. Not every combination
#' can be represented; a single probability-scale mode can remain unresolved
#' when the backend-default spread cannot support it as a unique global mode.
#'
#' ## Beliefs about residual standard deviation
#'
#' For Gaussian identity models, `update_priors("sigma", belief(mean = 30000))`
#' describes the prior mean of sigma itself. It does not fix sigma at 30000 or
#' describe a mean absolute residual. Use `sd` for uncertainty in sigma and a
#' 50% probability statement for its median. Only `scale = "original"` is
#' supported for sigma.
#'
#' Candidates are exponential, half-normal, half-Cauchy and half-Student's t,
#' using the positive zero-centered families supported by rstanarm's `prior_aux`
#' in version 2.32.2. All have mode zero; a mode outside the tolerance of zero
#' is backend-incompatible. No positive mode is reinterpreted as a mean.
#' Half-Cauchy has no finite mean or SD. Half-t mean requires df > 1 and SD
#' requires df > 2. These limitations apply when evaluating beliefs as well
#' as when displaying moments. A half-family scale is the underlying Normal
#' or t scale, not the SD of sigma's prior.
#'
#' A single mean, SD, or informative probability belief keeps the default
#' exponential family when compatible and determines its rate. With multiple
#' beliefs, the scale acceptance bands are intersected. The half-t search fits
#' scale and df separately from the coefficient t policy (fixed df = 7).
#' Its deterministic numerical search uses df above 0, 1, or 2 according to
#' requested moments, with an offset of at least `exp(-14)` and df at most
#' `1e6`. Failure within that finite search is not a proof of infeasibility.
#' All original belief tolerances must pass before a candidate is selected.
#'
#' Selection prefers a compatible exponential, then one-parameter half-normal
#' and half-Cauchy candidates, then the two-parameter half-t. This tie-breaker
#' does not establish uniqueness across families or within acceptance bands.
#' Uninformative beliefs such as `mode = 0` retain the disclosed default scale.
#' Diagnostics retain the original backend default. Expected-outcome
#' reconciliation never changes an elicited or default sigma prior.
#'
#' ## Beliefs about expected outcomes
#'
#' Supply the same `belief()` records to [update_outcome()] to describe a
#' conditional expected response at a chosen predictor configuration. This
#' currently requires a Gaussian identity model with an intercept, numeric
#' predictors, and independent Normal coefficient priors. Only original-scale
#' mean, mode, and probability beliefs are supported there, not SD beliefs.
#' These statements concern the expected response, not an individual noisy
#' observation. Reconciliation can change coefficient priors and compromise
#' earlier coefficient beliefs; their original statements remain available.
#' See [update_outcome()] for the reconciliation priorities and limitations.
#'
#' @return A `prior_belief` record that can be supplied to [update_priors()],
#'   [`update()`][update.flexbayes_priors], or [update_outcome()].
#' @seealso [identify_priors()], [plot.flexbayes_priors()], [check_prior()]
#' @examples
#' # These are separate statements, not one combined specification.
#' belief(mean = 100)
#' belief(mode = 100)
#' belief(sd = 20)
#' belief(value < 150, probability = .95)
#' belief(value > 200, probability = .10)
#'
#' # A median of 100 for a continuous distribution; no median argument needed.
#' belief(value < 100, probability = .50)
#'
#' # Thresholds can also use variables or ordinary numeric expressions.
#' upper = 150
#' belief(value <= upper, probability = .95)
#' belief(value > qt(.10, df = 7), probability = .90)
#'
#' # Each alternative scale has model and parameter restrictions (see Details).
#' belief(mode = .3, scale = "standardized")
#' belief(value < .8, probability = .95, scale = "probability")
#' belief(mean = 3, scale = "response")
#' belief(value < 2, probability = .95, scale = "odds_ratio")
#' belief(mode = 1.5, scale = "rate_ratio")
#'
#' # Model examples use the verified rstanarm 2.32.2 default-prior adapter.
#' if (requireNamespace("rstanarm", quietly = TRUE) &&
#'     utils::packageVersion("rstanarm") == "2.32.2") {
#'   d = data.frame(x = 1:6, y = c(200, 320, 470, 630, 760, 920))
#'   priors = identify_priors(y ~ x, data = d)
#'   priors = priors +
#'     update_priors(
#'       "x",
#'       belief(mode = 150),
#'       belief(value < 100, probability = .10)
#'     )
#'   print(summary(priors))
#'
#'   # Describe the expected response at x = 4, excluding observation noise.
#'   outcome_priors = identify_priors(y ~ x, data = d) +
#'     update_outcome(
#'       at = list(x = 4),
#'       belief(mode = 650),
#'       belief(value < 550, probability = .05)
#'     )
#'   print(summary(outcome_priors))
#' }
#' @export
belief = function(expression, probability = NULL, mean = NULL, mode = NULL, sd = NULL, tolerance = 1e-6, scale = "original") {
  # Check the requested units before interpreting the belief.
  if (length(scale) != 1L ||
      !is.character(scale) ||
      is.na(scale) ||
      !scale %in% c("original", "standardized", "probability", "response", "odds_ratio", "rate_ratio")) {
    stop("scale must be original, standardized, probability, response, odds_ratio or rate_ratio.")
  }

  # Keep the expression as entered, and identify any supplied moments.
  expr = substitute(expression)
  targets = list(mean = mean, mode = mode, sd = sd)
  supplied = !vapply(targets, is.null, logical(1))

  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) stop("Positive finite tolerance required.")

  # Resolve only the threshold, leaving the comparison itself unevaluated.
  if (!missing(expression)) {
    if (any(supplied) ||
        !is.call(expr) ||
        length(expr) != 3L ||
        !is.symbol(expr[[1]]) ||
        !as.character(expr[[1]]) %in% c("<", ">", "<=", ">=") ||
        !identical(expr[[2]], as.name("value"))) {
      stop("Use value < constant or value > constant.")
    }
    threshold = prior_resolve_threshold(expr[[3]], parent.frame())
    if (!is.numeric(probability) || length(probability) != 1 || !is.finite(probability) || probability < 0 || probability > 1) stop("Probability must be between zero and one.")
    out = list(
      type = "probability",
      target = probability,
      threshold = threshold,
      operator = as.character(expr[[1]]),
      expression = expr,
      label = paste(deparse(expr), "with probability", probability)
    )
  } else {
    if (sum(supplied) != 1L || !is.null(probability)) stop("Supply exactly one mean, mode, or sd.")
    type = names(targets)[supplied]
    target = targets[[type]]
    if (!is.numeric(target) || length(target) != 1 || !is.finite(target) || (type == "sd" && target <= 0)) stop("Invalid belief target.")
    out = list(type = type, target = target, label = paste(type, "=", target))
  }

  # Retain the entered units alongside the fitting tolerance.
  out$scale = scale
  if (scale != "original") out$label = paste0(out$label, " [", scale, "]")
  out$tolerance = tolerance
  class(out) = "prior_belief"
  out
}

prior_resolve_threshold = function(expression, environment) {
  threshold = tryCatch(eval(expression, envir = environment), error = function(error) {
    stop("Could not evaluate threshold: ", conditionMessage(error), call. = FALSE)
  })
  if (!is.numeric(threshold) || is.complex(threshold) || length(threshold) != 1L ||
      !is.null(dim(threshold)) || !is.finite(threshold)) {
    stop("Threshold must evaluate to a single finite numeric value.", call. = FALSE)
  }
  threshold
}

prior_fit_location_scale = function(records, reference, family = "normal") {
  # Repeated statements about the same quantity need overlapping tolerances.
  for (type in c("mean", "mode", "sd")) {
    same = Filter(function(r) r$type == type, records)
    if (length(same) > 1L) {
      lower = max(vapply(same, function(r) r$target - r$tolerance, numeric(1)))
      upper = min(vapply(same, function(r) r$target + r$tolerance, numeric(1)))
      if (lower > upper) return(list(valid = FALSE, message = paste("Inconsistent", type, "beliefs: target tolerances do not overlap.")))
    }
  }

  # Express each belief as a constraint on location and scale.
  entry = prior_distribution_registry()[[family]]
  unit = if (family == "normal") list(distribution = "normal", mean = 0, sd = 1) else list(distribution = "student_t", location = 0, scale = 1, df = prior_student_t_policy()$df)
  A = matrix(numeric(), ncol = 2)
  b = numeric()
  for (r in records) {
    row = switch(
      r$type,
      mean = c(1, 0),
      mode = c(1, 0),
      sd = c(0, entry$sd_factor),
      probability = c(1, prior_distribution_operation(unit, "quantile", if (r$operator %in% c("<", "<=")) r$target else 1 - r$target))
    )
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
  if (any(!is.finite(solution)) || solution[2] <= 0) {
    return(list(
      valid = FALSE,
      reference_used = rank < 2,
      message = "No location-scale prior with positive scale satisfies these beliefs."
    ))
  }
  distribution = if (family == "normal") list(distribution = family, mean = solution[1], sd = solution[2]) else list(
    distribution = family,
    location = solution[1],
    scale = solution[2],
    df = prior_student_t_policy()$df
  )

  # Check every original belief before accepting the fitted distribution.
  achieved = vapply(records, prior_distribution_target, numeric(1), prior = distribution)
  residuals = (achieved - vapply(records, `[[`, numeric(1), "target")) / vapply(records, `[[`, numeric(1), "tolerance")
  valid = all(abs(residuals) <= 1)
  list(
    valid = valid,
    reference_used = rank < 2,
    distribution = distribution,
    achieved = achieved,
    residuals = residuals,
    score = sum(residuals^2),
    rank = if (valid) 1L else NA_integer_,
    message = if (!valid) "Beliefs cannot be satisfied by this distribution family within tolerance." else NULL
  )
}

prior_distribution_target = function(record, prior) {
  if (prior_outcome_scale(record$scale)) return(prior_outcome_target(record, prior))
  if (record$type != "probability") return(prior_distribution_operation(prior, record$type))
  entry = prior_distribution_registry()[[prior$distribution]]
  if (record$operator %in% c(">", ">=") && is.function(entry$survival)) {
    return(prior_distribution_operation(prior, "survival", record$threshold))
  }
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
  # Validate additions and removals before changing the returned specification.
  if (length(parameter) != 1 || !parameter %in% names(object$beliefs)) stop("Unknown parameter.")
  prior_validate_standardized_beliefs(object)
  additions = list(...)
  if (any(!vapply(additions, inherits, logical(1), "prior_belief"))) stop("All additions must be belief objects.")
  if (parameter == "sigma" && prior_has_sigma(object) &&
      any(vapply(additions, function(r) r$scale != "original", logical(1)))) {
    stop("Sigma beliefs require scale = 'original': they describe the residual standard deviation itself.")
  }
  records = object$beliefs[[parameter]]
  if (!is.numeric(remove) ||
      anyNA(remove) ||
      any(!is.finite(remove)) ||
      any(remove != floor(remove)) ||
      any(remove < 1 |
      remove > length(records)) ||
      anyDuplicated(remove)) {
    stop("Invalid removal numbers.")
  }
  if (length(remove)) records = records[-remove]

  # Convert fitting units while keeping the original belief metadata.
  for (r in additions) {
    r = prior_convert_belief(r, object, parameter)
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

  # Mathematical contradictions take precedence over family fitting.
  consistency = if (parameter == "sigma" && prior_has_sigma(object))
    prior_scale_consistency(records) else prior_check_probability_constraints(records)
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
