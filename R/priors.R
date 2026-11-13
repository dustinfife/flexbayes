#' Identify Gaussian regression priors
#'
#' Defaults are resolved for rstanarm 2.32.2, dense design and QR = FALSE.
#' The intercept is at the column means of the model matrix, including
#' contrast and interaction columns. Defaults depend on mean(y) and sd(y).
#' No posterior fitting is performed.
#' @param formula Model formula.
#' @param data Data frame.
#' @param family Gaussian identity family.
#' @param package Backend name.
#' @return A flexbayes_priors object.
#' @export
identify_priors = function(formula, data, family = stats::gaussian(), package = "rstanarm") {
  if (!identical(package, "rstanarm")) stop("Only rstanarm is supported.")
  if (!requireNamespace("rstanarm", quietly = TRUE)) stop("Install rstanarm to resolve defaults.")
  version = as.character(utils::packageVersion("rstanarm"))
  if (version != "2.32.2") stop("Default resolution is verified only for rstanarm 2.32.2; installed: ", version)
  if (!is.list(family) || family$family != "gaussian" || family$link != "identity") stop("Only Gaussian identity models are supported.")
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("Supply a two-sided formula.")
  if (any(all.names(formula) %in% c("|", "||"))) stop("Mixed models are not supported.")
  frame = stats::model.frame(formula, data, na.action = stats::na.fail)
  terms = attr(frame, "terms")
  if (length(attr(terms, "offset"))) stop("Offsets are not supported yet.")
  X = stats::model.matrix(terms, frame)
  y = stats::model.response(frame)
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || nrow(X) < 2L || any(!is.finite(X))) stop("Finite numeric response and design required.")
  if (!ncol(X) || qr(X)$rank != ncol(X)) stop("A nonempty full-rank design is required.")
  sy = stats::sd(y)
  if (!is.finite(sy) || sy <= 0) stop("Backend default scaling requires positive sd(y).")
  intercept = attr(terms, "intercept") == 1L
  slopes = setdiff(colnames(X), "(Intercept)")
  reference = if (intercept) colMeans(X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  defaults = list()
  for (name in colnames(X)) {
    predictor_scale = if (length(unique(X[, name])) == 1L) 1 else stats::sd(X[, name])
    scale = if (name == "(Intercept)") 2.5 * sy else max(1e-12, 2.5 * sy / predictor_scale)
    scale = min(.Machine$double.xmax, scale)
    defaults[[name]] = list(distribution = "normal", mean = if (name == "(Intercept)") mean(y) else 0, sd = scale)
  }
  if ("sigma" %in% colnames(X)) stop("Predictor name sigma conflicts with the residual scale parameter.")
  defaults$sigma = list(distribution = "exponential", rate = 1 / sy)
  parameters = data.frame(parameter = names(defaults), role = c(ifelse(colnames(X) == "(Intercept)", "centered intercept", "coefficient"), "residual scale"), term = c(attr(X, "assign"), NA_integer_), lower = c(rep(-Inf, ncol(X)), 0), upper = Inf)
  out = list(model = list(formula = formula, data = data, family = family, frame = frame, terms = terms, X = X, contrasts = attr(X, "contrasts"), reference = reference), parameters = parameters, defaults = defaults, selected = defaults, beliefs = stats::setNames(rep(list(list()), length(defaults)), names(defaults)), candidates = stats::setNames(rep(list(list()), length(defaults)), names(defaults)), next_id = 1L, backend = list(package = package, version = version, compatible = TRUE, QR = FALSE, sparse = FALSE, outcome_mean = mean(y), outcome_sd = sy))
  out$model$standardization = prior_standardization_metadata(out$model)
  class(out) = "flexbayes_priors"
  out
}

prior_format = function(prior) {
  if (is.null(prior)) return("None selected")
  switch(prior$distribution,
    normal = sprintf("Normal(%.2f, %.2f)", prior$mean, prior$sd),
    student_t = sprintf("t(df=%g, %.2f, %.2f)", prior$df, prior$location, prior$scale),
    exponential = sprintf("Exponential(%.3g)", prior$rate),
    "Unsupported")
}

prior_console_line = function(text) {
  cat(paste(strwrap(text, width = max(40L, getOption("width", 80L))), collapse = "\n"), "\n", sep = "")
}

#' @export
print.flexbayes_priors = function(x, ...) {
  cat("Priors (", x$backend$package, " ", x$backend$version, ")\n", sep = "")
  source = vapply(names(x$selected), function(name) {
    if (identical(x$candidates[[name]]$status, "contradictory")) return("Contradictory")
    if (identical(x$candidates[[name]]$status, "backend_incompatible")) return("Incompatible")
    if (is.null(x$selected[[name]])) return("Unresolved")
    if (!length(x$beliefs[[name]])) return("Backend default")
    if (isTRUE(x$candidates[[name]]$reference_used)) return("Elicited*")
    "Elicited"
  }, character(1))
  table = data.frame(Parameter = names(x$selected),
    Prior = vapply(x$selected, prior_format, character(1)),
    Source = unname(source), Beliefs = lengths(x$beliefs), row.names = NULL)
  print(table, row.names = FALSE, right = FALSE)
  if (any(source == "Elicited*")) prior_console_line("* Backend default used to resolve underdetermined beliefs.")
  invisible(x)
}

#' Summarize prior specifications
#' @param object A prior specification.
#' @param diagnostics Include constraint targets, achieved values, standardized
#'   residuals, and candidate fit scores.
#' @param ... Reserved.
#' @return The original object, invisibly, after displaying its summary.
#' @export
summary.flexbayes_priors = function(object, diagnostics = FALSE, ...) {
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) stop("diagnostics must be TRUE or FALSE.")
  print(object)
  for (name in names(object$beliefs)) {
    records = object$beliefs[[name]]
    if (!length(records)) next
    cat("\nBeliefs for ", name, ":\n", sep = "")
    for (i in seq_along(records)) prior_console_line(paste0("  ", i, ". ", records[[i]]$label))
    if (diagnostics) {
      for (i in seq_along(records)) {
        r = records[[i]]
        if (!identical(r$scale, "standardized")) next
        equivalent = if (r$type == "probability") paste0("value ", r$operator, " ", format(r$threshold, digits = 6), ", probability = ", r$target) else paste(r$type, "=", format(r$target, digits = 6))
        prior_console_line(paste0("Belief ", i, ": SD(y)/SD(x) = ", format(r$standardization_factor, digits = 6), "; original-scale equivalent: ", equivalent))
      }
    }
    candidate = object$candidates[[name]]
    if (!is.null(candidate$message)) prior_console_line(paste0("Warning (", name, "): ", candidate$message))
    if (diagnostics && identical(candidate$status, "contradictory")) {
      cat("  Candidate fitting skipped: mathematical contradiction.\n")
    } else if (diagnostics) {
      cat("  ", if (identical(candidate$family, "student_t")) "Student's t" else "Normal", " candidate: ", if (isTRUE(candidate$valid)) "accepted" else "rejected", "\n", sep = "")
      if (length(candidate$alternatives)) {
        for (alternative in candidate$alternatives) {
          prior_console_line(paste0("  ", alternative$family, ": ", if (isTRUE(alternative$valid)) paste0("acceptable, mathematical rank ", alternative$rank) else paste0("rejected; ", alternative$message)))
        }
        if (!is.null(candidate$mathematical)) prior_console_line(paste("Mathematical preference:", prior_format(candidate$mathematical)))
      }
      if (!is.null(candidate$score)) cat("  Score: ", format(candidate$score, digits = 4), "\n", sep = "")
      if (!is.null(candidate$achieved)) {
        details = data.frame(Belief = seq_along(records),
          Target = vapply(records, `[[`, numeric(1), "target"),
          Achieved = candidate$achieved,
          Tolerance = vapply(records, `[[`, numeric(1), "tolerance"),
          Residual = candidate$residuals)
        print(details, row.names = FALSE, digits = 4)
        prior_console_line("Residuals are divided by each belief's tolerance; acceptable absolute values are <= 1.")
      } else prior_console_line("Constraint diagnostics unavailable: candidate rejected before evaluation.")
    }
  }
  cat("\n")
  prior_console_line("Backend defaults use observed mean(y) and sd(y).")
  if ("(Intercept)" %in% names(object$selected)) {
    prior_console_line("Centered intercept reference: model-matrix column means (including contrast and interaction columns).")
    reference = object$model$reference
    if (length(reference)) print(data.frame(Column = names(reference), Reference = unname(reference)), row.names = FALSE, digits = 4) else cat("Intercept-only model; no predictor columns.\n")
  } else cat("No intercept; predictors are not centered.\n")
  invisible(object)
}

#' Plot selected and candidate prior densities
#' @param x Prior specification.
#' @param parameter Optional parameter name. When supplied, compare acceptable candidates.
#' @param scale Plot in original or standardized slope units. Standardized plots require one supported parameter.
#' @param ... Reserved.
#' @return A ggplot object.
#' @export
plot.flexbayes_priors = function(x, parameter = NULL, scale = "original", ...) {
  if (length(scale) != 1L || !is.character(scale) || is.na(scale) || !scale %in% c("original", "standardized")) stop("scale must be original or standardized.")
  factor = 1
  if (scale == "standardized") {
    if (is.null(parameter) || length(parameter) != 1L) stop("Select one parameter for a standardized plot.")
    factor = prior_standardization_factor(x, parameter)
  }
  names = if (is.null(parameter)) names(x$selected) else parameter
  if (!length(names) || any(!names %in% names(x$selected))) stop("Unknown parameter.")
  rows = list()
  for (name in names) {
    distributions = list(Selected = x$selected[[name]])
    if (!is.null(parameter)) {
      alternatives = Filter(function(c) isTRUE(c$valid), x$candidates[[name]]$alternatives)
      distributions = lapply(alternatives, `[[`, "distribution")
      if (!length(distributions)) distributions = list(Selected = x$selected[[name]])
    }
    distributions = Filter(Negate(is.null), distributions)
    if (!length(distributions)) stop("No selected prior for ", name)
    for (label in names(distributions)) {
      p = distributions[[label]]
      grid = prior_distribution_operation(p, "quantile", seq(.001, .999, length.out = 501))
      density = prior_distribution_operation(p, "density", grid)
      selected = identical(p, x$selected[[name]])
      rows[[length(rows) + 1L]] = data.frame(value = grid / factor, density = density * factor, parameter = name,
        Candidate = paste0(label, if (selected) " (selected)" else " (mathematical only)"))
    }
  }
  ggplot2::ggplot(do.call(rbind, rows), ggplot2::aes(x = value, y = density, colour = Candidate)) + ggplot2::geom_line() + ggplot2::facet_wrap(~parameter, scales = "free") + ggplot2::labs(x = if (scale == "standardized") "Standardized coefficient" else "Parameter value", y = "Prior density")
}
