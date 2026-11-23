#' Switch the modeling package of an existing prior specification
#'
#' Preserve the original model, universal proper defaults, and all parameter and
#' expected-outcome beliefs while reevaluating priors for another backend. Belief
#' values, scales, tolerances, displayed order and stable internal IDs are kept;
#' threshold expressions are not evaluated again. No model is fitted or sampled.
#'
#' All parameter candidates are refitted before model-wide compatibility and
#' expected-outcome reconciliation are recomputed. Switching to brms enables its
#' additional sigma families and coefficient-specific priors. Switching to
#' rstanarm enforces its common slope-family and residual-prior restrictions.
#'
#' Compatible previous priors remain available as fallback references. A prior
#' unsupported by the destination is replaced by its original universal default
#' as a fallback reference. Mixed slope families are normalized to Normal
#' references before rstanarm selection, retaining existing Normal references.
#' These references do not override beliefs: acceptable candidates are selected
#' by the existing rules. Unsuccessful elicitation retains all beliefs, warns,
#' and still requires allow_fallback = TRUE for export, fitting or simulation.
#' The original object is never modified. Changed data or stale reconciliation
#' metadata must be resolved before switching. Calling with the current backend
#' also refits candidates and retries reconciliation.
#'
#' Belief updates automatically try brms only when rstanarm restrictions block
#' selection. The switch is accepted only after the entire specification, including
#' outcome reconciliation, passes brms export validation. Contradictions and
#' unverified numerical failures do not trigger migration. No automatic switch
#' back is made, even after removing beliefs. An explicit switch_package() call
#' evaluates exactly the requested package, without automatic reversal.
#' Automatic migrations are recorded in package_switches and displayed by
#' summary(..., diagnostics = TRUE).
#' @param object A flexbayes_priors object.
#' @param package Destination backend, either "rstanarm" or "brms".
#' @return A revised flexbayes_priors object with refreshed backend metadata,
#'   candidate rankings, compatible selections and outcome diagnostics.
#' @seealso [identify_priors()], [update_priors()], [as_rstanarm()], [as_brms()]
#' @examples
#' \dontrun{
#' priors = identify_priors(sale_price ~ sqft + age + bedrooms, homes_2024) +
#'   update_priors("sqft", belief(mode = 212))
#' # An optional explicit choice; compatible beliefs do not require a switch.
#' priors = switch_package(priors, package = "brms")
#' summary(priors)
#' check = check_prior(priors)
#' }
#' @export
switch_package = function(object, package) {
  if (!inherits(object, "flexbayes_priors")) stop("Supply a flexbayes_priors object.")
  version = prior_backend_version(package)
  prior_switch_validate(object)

  # Rebind only backend metadata; the estimation sample and reference priors stay fixed.
  object$backend$package = package
  object$backend$version = version
  object$backend$coefficients = NULL
  if (package == "brms") object$backend$coefficients = prior_brms_design(object)$coefficients
  prior_validate_response(stats::model.response(object$model$frame), object$model$family, package)
  object$retention_model = prior_retention_model(object)
  object = prior_switch_references(object)
  object$joint = NULL
  object$fallbacks = list()

  # Fit the complete belief collection before selecting any destination priors.
  for (parameter in names(object$beliefs)) object$candidates[[parameter]] = prior_parameter_candidates(object, parameter)
  object = prior_select_compatible(object)
  issues = prior_fallback_issues(object)
  if (length(issues)) warning("Switched to ", package, " with unsuccessful elicitation: ",
    paste(issues, collapse = ", "), ". All beliefs are preserved; inspect summary(). ",
    "Target-compatible fallbacks require explicit allow_fallback = TRUE for use.", call. = FALSE)
  object
}

prior_backend_version = function(package) {
  if (!is.character(package) || length(package) != 1L || !package %in% c("rstanarm", "brms")) stop("Choose package = 'rstanarm' or 'brms'.")
  if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " to identify or switch priors.")
  version = as.character(utils::packageVersion(package))
  if (package == "rstanarm" && version != "2.32.2") stop("Default resolution is verified only for rstanarm 2.32.2; installed: ", version)
  version
}

prior_switch_supported = function(parameter, object) {
  prior = object$selected[[parameter]]
  families = if (parameter == "sigma" && prior_has_sigma(object))
    prior_scale_families(object$backend$package) else c("normal", "student_t")
  !is.null(prior) && prior$distribution %in% families
}

prior_switch_validate = function(object) {
  if (!isTRUE(object$backend$compatible) || !prior_available_compatible(object)) stop("Current selected priors are not a coherent backend-compatible collection; identify priors again.")
  if (!identical(object$retention_model, prior_retention_model(object))) stop("Model parameterization or backend metadata changed; identify priors again.")
  prior_validate_family(object$model$family)
  prior_validate_standardization_data(object)
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  reference = if ("(Intercept)" %in% colnames(object$model$X)) colMeans(object$model$X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  if (!isTRUE(all.equal(reference, object$model$reference))) stop("Centering reference differs from the recorded model.")
  for (parameter in prior_required_parameters(object)) {
    if (!prior_switch_supported(parameter, object)) stop("Current prior for ", parameter, " is unsupported by its recorded backend; identify priors again.")
  }
  if (length(prior_fallback_issues(object))) prior_validate_retained(object)
  if (!isTRUE(object$joint$retained)) prior_joint_validate(object)
  invisible(object)
}

prior_switch_references = function(object) {
  reset = names(object$selected)[!vapply(names(object$selected), prior_switch_supported, logical(1), object = object)]
  slopes = setdiff(colnames(object$model$X), "(Intercept)")
  if (object$backend$package == "rstanarm") {
    families = vapply(object$selected[slopes], `[[`, character(1), "distribution")
    if (length(unique(families)) > 1L) reset = union(reset, slopes[families != "normal"])
  }
  for (parameter in reset) object$selected[parameter] = list(object$defaults[[parameter]])
  object$last_valid = object$selected
  object$last_valid_source[reset] = "backend_default"
  object$backend$compatible = prior_available_compatible(object)
  object
}

# Attempt migration once, only after selection has identified a package restriction.
# A trial is ordinary R data; failed trials leave the original fallback intact.
prior_maybe_switch_package = function(object) {
  if (object$backend$package != "rstanarm") return(object)
  statuses = vapply(object$candidates, function(candidate) if (is.null(candidate$status)) "default" else candidate$status, character(1))
  blocked = names(statuses)[statuses == "backend_incompatible"]
  if (!length(blocked) || any(statuses %in% c("contradictory", "no_supported_fit")) ||
      identical(object$joint$status, "incompatible")) return(object)
  if (!prior_package_available("brms")) {
    message("rstanarm cannot accommodate the requested prior for ", paste(blocked, collapse = ", "),
      ". Install brms to evaluate an alternative; retained priors still require allow_fallback = TRUE.")
    return(object)
  }

  # Refitting, reconciliation and faithful export must all pass before accepting.
  destination = tryCatch(suppressWarnings({
    candidate = switch_package(object, package = "brms")
    if (length(prior_fallback_issues(candidate))) NULL else {
      as_brms(candidate)
      candidate
    }
  }), error = function(error) NULL)
  if (is.null(destination)) return(object)
  reason = paste(blocked, collapse = ", ")
  destination$package_switches = c(object$package_switches,
    list(list(from = "rstanarm", to = "brms", parameters = blocked,
      reason = "rstanarm prior-family restrictions")))
  message("rstanarm cannot accommodate the requested prior for ", reason,
    ". Switching to brms after verifying the complete prior specification.")
  destination
}

prior_package_available = function(package) requireNamespace(package, quietly = TRUE)

prior_package_switch_summary = function(object) {
  for (event in object$package_switches) prior_console_line(paste0("Automatic package switch: ",
    event$from, " -> ", event$to, "; ", paste(event$parameters, collapse = ", "),
    " (", event$reason, ")."))
  invisible(NULL)
}
