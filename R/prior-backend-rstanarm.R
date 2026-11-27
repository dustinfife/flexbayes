#' Export a prior specification to rstanarm
#' @param object A flexbayes_priors object.
#' @param allow_fallback Use retained priors despite unsuccessful elicitation?
#'   Defaults to FALSE. TRUE explicitly approves their use and issues a warning;
#'   current beliefs remain recorded and are not treated as satisfied. Model-data
#'   validation and backend restrictions still apply.
#' @return A named argument list for rstanarm::stan_glm, containing formula,
#'   data, family, contrasts, na.action, prior, prior_intercept, prior_aux,
#'   QR and sparse. prior_aux is omitted for binomial/Poisson models, which
#'   have no auxiliary parameter. All resolved priors disable autoscaling. The intercept
#'   prior applies at the original model-matrix column means. Gaussian sigma
#'   priors export as exponential, zero-centered normal, or zero-centered t
#'   (df = 1 for half-Cauchy), with the backend's positive restriction and
#'   autoscaling disabled.
#' @export
as_rstanarm = function(object, allow_fallback = FALSE) {
  prior_rstanarm_validate(object, allow_fallback)

  # Export resolved priors without asking the backend to rescale them.
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  arguments = list(
    formula = object$model$formula,
    data = object$model$data,
    family = object$model$family,
    contrasts = object$model$contrasts,
    na.action = stats::na.fail,
    prior = prior_rstanarm_coefficients(object$selected[slopes]),
    prior_intercept = if ("(Intercept)" %in% colnames(object$model$X)) prior_rstanarm_coefficients(object$selected["(Intercept)"]) else NULL,
    prior_aux = if (prior_has_sigma(object)) prior_rstanarm_auxiliary(object$selected$sigma) else NULL,
    QR = FALSE,
    sparse = FALSE
  )

  # No auxiliary parameter exists for these GLMs. Omit the argument: an
  # explicit NULL is rejected by rstanarm prior-only mode even for GLMs.
  if (!prior_has_sigma(object)) arguments$prior_aux = NULL
  arguments
}

prior_rstanarm_validate = function(object, allow_fallback = FALSE) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  prior_fallback_permission(object, allow_fallback)
  if (object$backend$package != "rstanarm") stop("Only rstanarm export is supported.")
  if (!requireNamespace("rstanarm", quietly = TRUE)) stop("Install rstanarm to export priors.")
  version = as.character(utils::packageVersion("rstanarm"))
  if (version != "2.32.2" || object$backend$version != version) stop("Faithful export is verified only for matching rstanarm 2.32.2 specifications.")
  prior_validate_family(object$model$family)

  # Check compatibility across the complete collection of slope priors.
  blocked = names(Filter(function(c) identical(c$status, "backend_incompatible"), object$candidates))
  if (length(blocked) && !allow_fallback && !length(prior_fallback_issues(object))) stop("Incompatible parameters: ", paste(blocked, collapse = ", "),
    ". rstanarm requires compatible prior_aux distributions and one slope prior family across the model.")
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  selected = Filter(Negate(is.null), object$selected[slopes])
  families = unique(vapply(selected, `[[`, character(1), "distribution"))
  if (length(families) > 1L) stop("Incompatible coefficient families: ", paste(names(selected), collapse = ", "), ". rstanarm requires one slope prior family.")
  prior_predictive_validate(object, allow_fallback)
  for (name in names(object$beliefs)) {
    records = Filter(function(r) !prior_is_r_squared(r), object$beliefs[[name]])
    if (!length(records)) next
    if (allow_fallback && !is.null(object$fallbacks[[name]])) next
    achieved = vapply(records, prior_distribution_target, numeric(1), prior = object$selected[[name]])
    targets = vapply(records, `[[`, numeric(1), "target")
    tolerances = vapply(records, `[[`, numeric(1), "tolerance")
    if ((is.null(object$joint) ||
        name == "sigma") &&
        (any(!is.finite(achieved)) ||
        any(abs(achieved - targets) > tolerances))) {
      stop("Selected prior does not satisfy recorded beliefs for ", name, ".")
    }
  }
  if (prior_has_sigma(object) && !object$selected$sigma$distribution %in% prior_scale_families()) {
    stop("Unsupported residual-scale prior family for rstanarm prior_aux.")
  }
  if (!identical(object$backend$QR, FALSE) || !identical(object$backend$sparse, FALSE)) stop("Export requires QR = FALSE and sparse = FALSE.")

  # Reject changed data instead of exporting priors with a new interpretation.
  frame = stats::model.frame(object$model$formula, object$model$data, na.action = stats::na.fail)
  X = stats::model.matrix(attr(frame, "terms"), frame, contrasts.arg = object$model$contrasts)
  if (!isTRUE(all.equal(X, object$model$X)) ||
      !isTRUE(all.equal(stats::model.response(frame), stats::model.response(object$model$frame)))) {
    stop("Model data or formula evaluation changed since prior identification; identify priors again.")
  }
  slopes = setdiff(colnames(X), "(Intercept)")
  reference = if ("(Intercept)" %in% colnames(X)) colMeans(X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  if (!isTRUE(all.equal(reference, object$model$reference))) stop("Centering reference differs from the recorded model.")
  invisible(object)
}

prior_rstanarm_coefficients = function(priors) {
  if (!length(priors)) return(NULL)
  family = unique(vapply(priors, `[[`, character(1), "distribution"))
  if (length(family) != 1L) stop("Incompatible coefficient families: ", paste(names(priors), collapse = ", "), "; rstanarm requires one slope prior family.")
  if (family == "normal") {
    return(rstanarm::normal(
      location = vapply(priors, `[[`, numeric(1), "mean"),
      scale = vapply(priors, `[[`, numeric(1), "sd"), autoscale = FALSE))
  }
  if (family == "student_t") {
    return(rstanarm::student_t(
      df = vapply(priors, `[[`, numeric(1), "df"),
      location = vapply(priors, `[[`, numeric(1), "location"),
      scale = vapply(priors, `[[`, numeric(1), "scale"), autoscale = FALSE))
  }
  stop("Unsupported coefficient prior family: ", family)
}

prior_rstanarm_fit = function(arguments) {
  do.call(rstanarm::stan_glm, arguments)
}

# continuous.stan in rstanarm 2.32.2 restricts aux_unscaled to positive values.
# Zero location is essential: nonzero locations do not describe these halves.
prior_rstanarm_auxiliary = function(prior) {
  prior_validate_distribution(prior)
  switch(prior$distribution,
    exponential = rstanarm::exponential(rate = prior$rate, autoscale = FALSE),
    half_normal = rstanarm::normal(location = 0, scale = prior$scale, autoscale = FALSE),
    half_student_t = rstanarm::student_t(df = prior$df, location = 0, scale = prior$scale, autoscale = FALSE),
    half_cauchy = rstanarm::cauchy(location = 0, scale = prior$scale, autoscale = FALSE),
    stop("Unsupported residual-scale prior family for rstanarm prior_aux."))
}

# Effective defaults verified against stan_glm.fit in rstanarm 2.32.2.
prior_rstanarm_defaults = function(X, y, family) {
  gaussian = family$family == "gaussian"
  sy = if (gaussian) stats::sd(y) else 1
  defaults = list()

  # Mirror the verified backend defaults, including its scaling safeguards.
  for (name in colnames(X)) {
    sx = if (length(unique(X[, name])) == 1L) 1 else stats::sd(X[, name])
    scale = if (name == "(Intercept)") 2.5 * sy else max(1e-12, 2.5 * sy / sx)
    defaults[[name]] = list(distribution = "normal", mean = if (gaussian && name == "(Intercept)") mean(y) else 0,
      sd = min(.Machine$double.xmax, scale))
  }
  if (gaussian) defaults$sigma = list(distribution = "exponential", rate = 1 / sy)
  defaults
}
