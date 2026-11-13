#' Fit a model using its selected priors
#' @param object A prior specification.
#' @param ... Named backend sampling or algorithm arguments. Model and prior
#'   arguments cannot be overridden. See details.
#' @details Supported arguments are chains, iter, warmup, thin, seed, cores,
#' refresh, control, init, save_warmup, sample_file, diagnostic_file, verbose,
#' algorithm, adapt_delta, mean_PPD and prior_PD. Unknown and partially matched
#' names are rejected. prior_PD = TRUE explicitly requests prior-only sampling.
#' @return The backend's native fitted model (stanreg for rstanarm).
#' @export
fit = function(object, ...) UseMethod("fit")

#' @rdname fit
#' @export
fit.flexbayes_priors = function(object, ...) {
  arguments = list(...)
  allowed = c("chains", "iter", "warmup", "thin", "seed", "cores", "refresh",
    "control", "init", "save_warmup", "sample_file", "diagnostic_file", "verbose",
    "algorithm", "adapt_delta", "mean_PPD", "prior_PD")
  names = names(arguments)
  if (length(arguments) && (is.null(names) || anyNA(names) || any(!nzchar(names)) || anyDuplicated(names))) stop("Fitting arguments must have unique, complete names.")
  prohibited = setdiff(names, allowed)
  if (length(prohibited)) stop("Cannot override model, priors, or parameterization; unsupported fitting arguments: ", paste(prohibited, collapse = ", "))
  backend = object$backend$package
  if (!identical(backend, "rstanarm")) stop("Unsupported fitting backend: ", backend)
  prior_rstanarm_fit(c(as_rstanarm(object), arguments))
}
