#' Export a prior specification to brms
#'
#' Supplies explicit proper priors for every supported parameter. These starting
#' defaults are flexbayes priors derived from rstanarm 2.32.2, not native brms
#' defaults. Coefficient priors can use different families. The intercept prior
#' applies at the original model-matrix column means, including contrast and
#' interaction columns. Binary responses export as Bernoulli logit models.
#' @param object A flexbayes_priors object currently using brms, whether selected
#'   initially, explicitly with switch_package(), or automatically after an update.
#' @param allow_fallback Explicitly approve retained priors after unsuccessful
#'   elicitation, with a warning. Defaults to FALSE; beliefs remain recorded.
#' @return A named argument list for brms::brm containing formula, data, family,
#'   and prior. The data preserve the contrasts used during identification.
#' @seealso [identify_priors()], [as_rstanarm()], [fit()]
#' @export
as_brms = function(object, allow_fallback = FALSE) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  if (object$backend$package != "brms") stop("Only brms specifications can be exported with as_brms().")
  prior_predictive_validate(object, allow_fallback)
  prior_brms_validate_beliefs(object, allow_fallback)
  design = prior_brms_design(object)
  if (!identical(design$coefficients, object$backend$coefficients)) stop("brms coefficient mapping changed; identify priors again.")

  # Set each coefficient explicitly; never inherit brms's flat slope defaults.
  priors = list()
  for (name in names(object$selected)) {
    class = if (name == "(Intercept)") "Intercept" else if (name == "sigma" && prior_has_sigma(object)) "sigma" else "b"
    coefficient = if (class == "b") unname(design$coefficients[name]) else ""
    priors[[name]] = brms::set_prior(prior_brms_distribution(object$selected[[name]]), class = class, coef = coefficient)
  }
  list(formula = object$model$formula, data = design$data,
    family = prior_brms_family(object), prior = do.call(c, unname(priors)))
}

prior_brms_family = function(object) {
  if (object$model$family$family == "binomial") return(brms::bernoulli(link = "logit"))
  object$model$family
}

prior_brms_data = function(object) {
  data = object$model$data
  for (name in names(object$model$contrasts)) {
    if (!name %in% names(data) || !is.factor(data[[name]])) next
    stats::contrasts(data[[name]]) = object$model$contrasts[[name]]
  }
  data
}

prior_brms_validate_data = function(object) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("Install brms to use this specification.")
  if (as.character(utils::packageVersion("brms")) != object$backend$version) stop("brms version changed; identify priors again.")
  prior_validate_family(object$model$family)
  prior_validate_standardization_data(object)
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  reference = if ("(Intercept)" %in% colnames(object$model$X)) colMeans(object$model$X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  if (!isTRUE(all.equal(reference, object$model$reference))) stop("Centering reference differs from the recorded model.")
  if (!identical(object$backend$QR, FALSE) || !identical(object$backend$sparse, FALSE)) stop("brms export requires the original dense design without QR.")
  invisible(object)
}

prior_brms_design = function(object) {
  data = prior_brms_data(object)
  family = prior_brms_family(object)
  available = brms::get_prior(object$model$formula, data = data, family = family)
  standata = brms::make_standata(object$model$formula, data = data, family = family)
  X = standata$X
  original = object$model$X
  if (is.null(X) || !identical(dim(X), dim(original)) ||
      !isTRUE(all.equal(unname(X), unname(original), check.attributes = FALSE))) {
    stop("brms model matrix differs from the identified design; this formula or parameterization is unsupported.")
  }
  intercept = "(Intercept)" %in% colnames(original)
  if (intercept && (!"Intercept" %in% available$class || colnames(X)[1] != "Intercept")) stop("brms centered-intercept parameterization could not be verified.")
  indices = which(colnames(original) != "(Intercept)")
  coefficients = stats::setNames(colnames(X)[indices], colnames(original)[indices])
  known = available$coef[available$class == "b" & nzchar(available$coef)]
  if (anyDuplicated(coefficients) || !setequal(unname(coefficients), known)) stop("brms coefficient names cannot be mapped uniquely to the identified model matrix.")
  list(data = data, coefficients = coefficients)
}

prior_brms_validate_beliefs = function(object, allow_fallback) {
  for (name in names(object$beliefs)) {
    if (allow_fallback && !is.null(object$fallbacks[[name]])) next
    if (!is.null(object$joint) && name != "sigma") next
    records = object$beliefs[[name]]
    if (!length(records)) next
    achieved = vapply(records, prior_distribution_target, numeric(1), prior = object$selected[[name]])
    target = vapply(records, `[[`, numeric(1), "target")
    tolerance = vapply(records, `[[`, numeric(1), "tolerance")
    if (any(!is.finite(achieved)) || any(abs(achieved - target) > tolerance)) stop("Selected prior does not satisfy recorded beliefs for ", name, ".")
  }
  invisible(object)
}

prior_brms_distribution = function(prior) {
  prior_validate_distribution(prior)
  values = switch(prior$distribution,
    normal = c("normal", prior$mean, prior$sd),
    student_t = c("student_t", prior$df, prior$location, prior$scale),
    exponential = c("exponential", prior$rate),
    half_normal = c("normal", 0, prior$scale),
    half_student_t = c("student_t", prior$df, 0, prior$scale),
    half_cauchy = c("cauchy", 0, prior$scale),
    gamma = c("gamma", prior$shape, prior$rate),
    lognormal = c("lognormal", prior$meanlog, prior$sdlog),
    stop("Unsupported brms prior distribution."))
  # Round-trip precision keeps export identical to simulation parameters.
  paste0(values[1], "(", paste(sprintf("%.17g", as.numeric(values[-1])), collapse = ", "), ")")
}

prior_brms_fit = function(arguments) do.call(brms::brm, arguments)
