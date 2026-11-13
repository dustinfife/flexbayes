prior_validate_family = function(family) {
  if (!is.list(family) || length(family$family) != 1L || length(family$link) != 1L ||
      !paste(family$family, family$link) %in% c("gaussian identity", "binomial logit", "poisson log")) stop("Supported families: Gaussian identity, binomial logit, and Poisson log only.")
  invisible(family)
}

prior_has_sigma = function(object) identical(object$model$family$family, "gaussian")
prior_required_parameters = function(object) c(colnames(object$model$X), if (prior_has_sigma(object)) "sigma")

prior_validate_response = function(y, family) {
  if (!is.numeric(y) || !is.null(dim(y)) || any(!is.finite(y))) stop("Response must be a finite numeric vector; binomial models require binary 0/1 outcomes.")
  if (family$family == "gaussian" && (!is.finite(stats::sd(y)) || stats::sd(y) <= 0)) stop("Backend default scaling requires positive sd(y).")
  if (family$family == "binomial" && any(!y %in% c(0, 1))) stop("Binomial logit models require binary 0/1 outcomes; grouped and proportion formats are not supported.")
  if (family$family == "poisson") {
    if (any(y < 0 | y != floor(y) | y > .Machine$integer.max)) stop("Poisson outcomes must be nonnegative integer counts within Stan's integer range.")
    if (length(unique(y)) == 1L && y[1] > 1) stop("rstanarm rejects constant responses greater than one; zero-count and all-zero responses are supported.")
  }
  invisible(y)
}

prior_simulate_response = function(linear, draws, family) {
  if (any(!is.finite(linear))) stop("Non-finite linear predictors; revise priors or predictor values.")
  if (family$family == "gaussian") return(list(expected = linear, simulated = prior_simulate_gaussian(linear, draws[, "sigma"])))
  expected = if (family$family == "binomial") stats::plogis(linear) else exp(linear)
  if (family$family == "poisson" && (any(!is.finite(expected)) || any(expected > 1e8))) stop("Poisson predicted rates are nonfinite or exceed the simulation safety limit of 1e8; revise priors or predictor values.")
  values = if (family$family == "binomial") stats::rbinom(length(expected), size = 1, prob = expected) else stats::rpois(length(expected), lambda = expected)
  list(expected = expected, simulated = matrix(values, nrow = nrow(linear), dimnames = dimnames(linear)))
}

prior_default_description = function(object) {
  if (prior_has_sigma(object)) return("Backend defaults use observed mean(y) and sd(y).")
  paste0("Priors use the ", if (object$model$family$family == "binomial") "log-odds" else "log-rate", " scale; default slope scales use predictor SDs, not outcome SD.")
}
