prior_distribution_registry = function() {
  list(normal = list(density = stats::dnorm, cdf = stats::pnorm,
    quantile = stats::qnorm, random = stats::rnorm, support = c(-Inf, Inf), sd_factor = 1),
    student_t = list(density = prior_dt, cdf = prior_pt, quantile = prior_qt,
      random = prior_rt, support = c(-Inf, Inf), sd_factor = sqrt(7 / 5)),
    exponential = list(density = stats::dexp, cdf = stats::pexp,
      quantile = stats::qexp, random = stats::rexp, support = c(0, Inf)))
}

prior_dt = function(x, location, scale, df) stats::dt((x - location) / scale, df) / scale
prior_pt = function(q, location, scale, df) stats::pt((q - location) / scale, df)
prior_qt = function(p, location, scale, df) location + scale * stats::qt(p, df)
prior_rt = function(n, location, scale, df) location + scale * stats::rt(n, df)

prior_validate_distribution = function(prior) {
  if (!is.list(prior) || length(prior$distribution) != 1L) stop("Invalid distribution specification.")
  required = switch(prior$distribution, normal = c("mean", "sd"),
    student_t = c("location", "scale", "df"), exponential = "rate", NULL)
  if (is.null(required)) stop("Unsupported distribution.")
  for (name in required) {
    value = prior[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) stop("Invalid distribution parameter: ", name)
  }
  positive = intersect(required, c("sd", "scale", "df", "rate"))
  if (any(vapply(prior[positive], function(v) v <= 0, logical(1)))) stop("Distribution scales, rates and df must be positive.")
  invisible(prior)
}

prior_distribution_operation = function(prior, operation, value = NULL) {
  prior_validate_distribution(prior)
  entry = prior_distribution_registry()[[prior$distribution]]
  if (operation == "support") return(entry$support)
  if (operation == "mode") return(switch(prior$distribution, normal = prior$mean, student_t = prior$location, exponential = 0))
  if (operation == "mean") return(switch(prior$distribution, normal = prior$mean, student_t = if (prior$df > 1) prior$location else NA_real_, exponential = 1 / prior$rate))
  if (operation == "sd") return(switch(prior$distribution, normal = prior$sd, student_t = if (prior$df > 2) prior$scale * sqrt(prior$df / (prior$df - 2)) else Inf, exponential = 1 / prior$rate))
  implementation = entry[[operation]]
  if (!is.function(implementation)) stop("Unsupported distribution operation.")
  arguments = prior[setdiff(names(prior), "distribution")]
  do.call(implementation, c(list(value), arguments))
}
