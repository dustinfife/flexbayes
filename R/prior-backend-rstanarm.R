#' Export a prior specification to rstanarm
#' @param object A flexbayes_priors object.
#' @return A named argument list for rstanarm::stan_glm, containing formula,
#'   data, family, contrasts, na.action, prior, prior_intercept, prior_aux,
#'   QR and sparse. prior_aux is omitted for binomial/Poisson models, which
#'   have no auxiliary parameter. All resolved priors disable autoscaling. The intercept
#'   prior applies at the original model-matrix column means.
#' @export
as_rstanarm = function(object) {
  prior_rstanarm_validate(object)
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  arguments = list(formula = object$model$formula, data = object$model$data,
    family = object$model$family, contrasts = object$model$contrasts,
    na.action = stats::na.fail,
    prior = prior_rstanarm_coefficients(object$selected[slopes]),
    prior_intercept = if ("(Intercept)" %in% colnames(object$model$X)) prior_rstanarm_coefficients(object$selected["(Intercept)"]) else NULL,
    prior_aux = if (prior_has_sigma(object)) rstanarm::exponential(rate = object$selected$sigma$rate, autoscale = FALSE) else NULL,
    QR = FALSE, sparse = FALSE)
  # No auxiliary parameter exists for these GLMs. Omit the argument: an
  # explicit NULL is rejected by rstanarm prior-only mode even for GLMs.
  if (!prior_has_sigma(object)) arguments$prior_aux = NULL
  arguments
}

prior_rstanarm_validate = function(object) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  if (object$backend$package != "rstanarm") stop("Only rstanarm export is supported.")
  if (!requireNamespace("rstanarm", quietly = TRUE)) stop("Install rstanarm to export priors.")
  version = as.character(utils::packageVersion("rstanarm"))
  if (version != "2.32.2" || object$backend$version != version) stop("Faithful export is verified only for matching rstanarm 2.32.2 specifications.")
  prior_validate_family(object$model$family)
  blocked = names(Filter(function(c) identical(c$status, "backend_incompatible"), object$candidates))
  if (length(blocked)) stop("Incompatible parameters: ", paste(blocked, collapse = ", "), ". rstanarm requires one slope prior family across the model.")
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  selected = Filter(Negate(is.null), object$selected[slopes])
  families = unique(vapply(selected, `[[`, character(1), "distribution"))
  if (length(families) > 1L) stop("Incompatible coefficient families: ", paste(names(selected), collapse = ", "), ". rstanarm requires one slope prior family.")
  prior_predictive_validate(object)
  for (name in names(object$beliefs)) {
    records = object$beliefs[[name]]
    if (!length(records)) next
    achieved = vapply(records, prior_distribution_target, numeric(1), prior = object$selected[[name]])
    targets = vapply(records, `[[`, numeric(1), "target")
    tolerances = vapply(records, `[[`, numeric(1), "tolerance")
    if (any(!is.finite(achieved)) || any(abs(achieved - targets) > tolerances)) stop("Selected prior does not satisfy recorded beliefs for ", name, ".")
  }
  if (prior_has_sigma(object) && object$selected$sigma$distribution != "exponential") stop("Only exponential residual-scale priors can currently be exported faithfully.")
  if (!identical(object$backend$QR, FALSE) || !identical(object$backend$sparse, FALSE)) stop("Export requires QR = FALSE and sparse = FALSE.")
  frame = stats::model.frame(object$model$formula, object$model$data, na.action = stats::na.fail)
  X = stats::model.matrix(attr(frame, "terms"), frame, contrasts.arg = object$model$contrasts)
  if (!isTRUE(all.equal(X, object$model$X)) ||
      !isTRUE(all.equal(stats::model.response(frame), stats::model.response(object$model$frame)))) stop("Model data or formula evaluation changed since prior identification; identify priors again.")
  slopes = setdiff(colnames(X), "(Intercept)")
  reference = if ("(Intercept)" %in% colnames(X)) colMeans(X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  if (!isTRUE(all.equal(reference, object$model$reference))) stop("Centering reference differs from the recorded model.")
  invisible(object)
}

prior_rstanarm_coefficients = function(priors) {
  if (!length(priors)) return(NULL)
  family = unique(vapply(priors, `[[`, character(1), "distribution"))
  if (length(family) != 1L) stop("Incompatible coefficient families: ", paste(names(priors), collapse = ", "), "; rstanarm requires one slope prior family.")
  if (family == "normal") return(rstanarm::normal(
    location = vapply(priors, `[[`, numeric(1), "mean"),
    scale = vapply(priors, `[[`, numeric(1), "sd"), autoscale = FALSE))
  if (family == "student_t") return(rstanarm::student_t(
    df = vapply(priors, `[[`, numeric(1), "df"),
    location = vapply(priors, `[[`, numeric(1), "location"),
    scale = vapply(priors, `[[`, numeric(1), "scale"), autoscale = FALSE))
  stop("Unsupported coefficient prior family: ", family)
}

prior_rstanarm_fit = function(arguments) {
  do.call(rstanarm::stan_glm, arguments)
}

# Effective defaults verified against stan_glm.fit in rstanarm 2.32.2.
prior_rstanarm_defaults = function(X, y, family) {
  gaussian = family$family == "gaussian"
  sy = if (gaussian) stats::sd(y) else 1
  defaults = list()
  for (name in colnames(X)) {
    sx = if (length(unique(X[, name])) == 1L) 1 else stats::sd(X[, name])
    scale = if (name == "(Intercept)") 2.5 * sy else max(1e-12, 2.5 * sy / sx)
    defaults[[name]] = list(distribution = "normal", mean = if (gaussian && name == "(Intercept)") mean(y) else 0,
      sd = min(.Machine$double.xmax, scale))
  }
  if (gaussian) defaults$sigma = list(distribution = "exponential", rate = 1 / sy)
  defaults
}
