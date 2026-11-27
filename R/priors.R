#' Identify regression priors
#'
#' Both backends start with the same proper flexbayes defaults, derived from
#' rstanarm 2.32.2 (dense design, QR = FALSE). Native brms defaults are deliberately
#' overridden; prior-origin metadata are stored separately from the fitting backend.
#' The intercept is at the column means of the model matrix, including
#' contrast and interaction columns. Gaussian defaults depend on mean(y) and
#' sd(y); GLM defaults use the link scale and predictor SDs.
#' No posterior fitting is performed.
#' @param formula Model formula.
#' @param data Data frame.
#' @param family Gaussian identity, binomial logit or Poisson log family.
#' @param package Initial modeling-package preference: rstanarm (default) or brms.
#'   Belief updates automatically switch from rstanarm to brms when package
#'   restrictions block selection and the complete specification can be verified
#'   and exported. brms must be installed; packages are never installed automatically.
#'   Removing beliefs does not switch back. Use switch_package() for explicit changes.
#' @details Unsuccessful updates preserve all beliefs and retain the last valid
#' prior, or the original backend default. Summary and plots show that retained
#' prior with a fallback label; they do not imply that current beliefs were met.
#' Failed outcome reconciliation retains the previous coherent coefficient-prior
#' collection. Export, fitting and simulation require explicit approval through
#' allow_fallback = TRUE before using these priors.
#' @return A flexbayes_priors object.
#' @seealso [switch_package()], [update_priors()]
#' @export
identify_priors = function(formula, data, family = stats::gaussian(), package = "rstanarm") {
  version = prior_backend_version(package)
  prior_validate_family(family)
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("Supply a two-sided formula.")
  if (any(all.names(formula) %in% c("|", "||"))) stop("Mixed models are not supported.")

  # Use one estimation sample for defaults, centering and later predictions.
  frame = stats::model.frame(formula, data, na.action = stats::na.fail)
  terms = attr(frame, "terms")
  if (length(attr(terms, "offset"))) stop("Offsets are not supported yet.")
  X = stats::model.matrix(terms, frame)
  y = stats::model.response(frame)
  if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) || nrow(X) < 2L || any(!is.finite(X))) stop("Finite numeric response and design required.")
  if (!ncol(X) || qr(X)$rank != ncol(X)) stop("A nonempty full-rank design is required.")
  prior_validate_response(y, family, package)
  sy = stats::sd(y)
  intercept = attr(terms, "intercept") == 1L
  slopes = setdiff(colnames(X), "(Intercept)")
  reference = if (intercept) colMeans(X[, slopes, drop = FALSE]) else stats::setNames(rep(0, length(slopes)), slopes)
  if (family$family == "gaussian" && "sigma" %in% colnames(X)) stop("Predictor name sigma conflicts with the residual scale parameter.")

  # Preserve backend defaults separately from any elicited selections.
  defaults = prior_rstanarm_defaults(X, y, family)
  parameters = data.frame(
    parameter = names(defaults),
    role = c(ifelse(colnames(X) == "(Intercept)", "centered intercept", "coefficient"), if (family$family == "gaussian") "residual scale"),
    term = c(attr(X, "assign"), if (family$family == "gaussian") NA_integer_),
    lower = c(rep(-Inf, ncol(X)), if (family$family == "gaussian") 0),
    upper = Inf
  )
  out = list(
    model = list(formula = formula, data = data, family = family, frame = frame, terms = terms, X = X, contrasts = attr(X, "contrasts"), reference = reference),
    parameters = parameters,
    default_origin = list(package = "rstanarm", version = "2.32.2", policy = "flexbayes proper defaults"),
    defaults = defaults,
    selected = defaults,
    beliefs = stats::setNames(rep(list(list()), length(defaults)), names(defaults)),
    candidates = stats::setNames(rep(list(list()), length(defaults)), names(defaults)),
    next_id = 1L,
    backend = list(package = package, version = version, compatible = TRUE, QR = FALSE, sparse = FALSE, outcome_mean = mean(y), outcome_sd = sy)
  )
  if (family$family != "gaussian") out$parameters$scale = if (family$family == "binomial") "log-odds" else "log-rate"

  if (family$family == "gaussian") out$model$predictor_covariance = prior_predictor_covariance(X)
  out$model$standardization = prior_standardization_metadata(out$model)
  class(out) = "flexbayes_priors"
  if (package == "brms") out$backend$coefficients = prior_brms_design(out)$coefficients
  prior_retention_initialize(out)
}

prior_format = function(prior) {
  if (is.null(prior)) return("None selected")
  switch(
    prior$distribution,
    normal = sprintf("Normal(%.2f, %.2f)", prior$mean, prior$sd),
    student_t = sprintf("t(df=%g, %.2f, %.2f)", prior$df, prior$location, prior$scale),
    exponential = sprintf("Exponential(%.3g)", prior$rate),
    half_normal = sprintf("Half-Normal(0, %.3g)", prior$scale),
    half_student_t = sprintf("Half-t(df=%.4g, 0, %.3g)", prior$df, prior$scale),
    half_cauchy = sprintf("Half-Cauchy(0, %.3g)", prior$scale),
    gamma = sprintf("Gamma(shape=%.4g, rate=%.4g)", prior$shape, prior$rate),
    lognormal = sprintf("Lognormal(%.4g, %.4g)", prior$meanlog, prior$sdlog),
    "Unsupported"
  )
}

prior_console_line = function(text) {
  cat(paste(strwrap(text, width = max(40L, getOption("width", 80L))), collapse = "\n"), "\n", sep = "")
}

#' @export
print.flexbayes_priors = function(x, ...) {
  cat("Priors (", x$backend$package, " ", x$backend$version, ")\n", sep = "")
  if (x$backend$package == "brms") prior_console_line("Defaults: flexbayes proper priors derived from rstanarm 2.32.2; native brms defaults are overridden.")
  indirect = prior_indirect_report(x)
  table = prior_report_table(x, indirect)
  source = table$Source
  print(table, row.names = FALSE, right = FALSE)
  if (any(source == "Elicited*")) prior_console_line(if (x$backend$package == "brms") "* Universal flexbayes default used to resolve underdetermined beliefs; see diagnostics for the retained reference moment." else "* Backend default used to resolve underdetermined beliefs.")
  prior_report_indirect_note(indirect)
  if (length(prior_fallback_issues(x))) prior_console_line("Fallback priors are available for inspection; current elicitation was unsuccessful. See summary() before using allow_fallback = TRUE.")
  if (!is.null(x$joint)) cat("Outcome reconciliation:", x$joint$status, ";", length(x$outcomes), "configuration(s).\n")
  invisible(x)
}

#' Summarize prior specifications
#'
#' Displays a compact prior table, original beliefs and relevant warnings. With
#' outcome beliefs, shows outcome satisfaction and one grouped comparison of
#' changed coefficient priors and their compromised beliefs. Before means the
#' original distribution or requested belief value; After means the adjusted
#' distribution or achieved value. Normal is abbreviated to N only in this table.
#' Unchanged priors and satisfied beliefs are omitted from the comparison.
#' Failed elicitation is labeled as a fallback, with its retained prior's origin
#' and the reason the current beliefs were not verified.
#' Unverified reconciliations are clearly flagged. Complete comparisons,
#' backend-centering references and optimization statistics require diagnostics=TRUE.
#' @details
#' The compact table in `print(priors)` distinguishes where each prior came
#' from. A changed prior informed by one outcome or difference specification
#' involving only that coefficient is labelled `Elicited via outcome` or
#' `Elicited via difference`. Shared constraints, several specifications, or
#' combined parameter and indirect beliefs use the general label `Reconciled`.
#' These labels describe the route, not a guarantee of unique identification.
#'
#' Plain belief counts refer to beliefs supplied for that parameter. An entry
#' such as `2 indirect` counts outcome/difference statements involving the
#' coefficient. Mixed counts such as `1 + 2` mean one parameter belief plus two
#' indirect statements. A shared statement can involve several coefficients:
#' row counts must not be added as if they were separate pieces of information.
#' The total counts each supplied statement once, including repeated statements.
#'
#' Participation does not imply a changed prior. Unchanged participants keep
#' their previous source; the grouped summary lists them separately. Changed
#' groups display the same source and counts as `print(priors)`. Unverified
#' indirect participation is labelled pending, and existing fallback labels
#' remain in place. R-squared beliefs supplied for sigma stay in sigma's
#' parameter-belief count, not the indirect counts for coefficients.
#' @param object A prior specification.
#' @param diagnostics Include constraint targets, achieved values, standardized
#'   residuals, and candidate fit scores. Outcome reconciliation additionally
#'   reports configurations, original and adjusted coefficient priors, relaxed
#'   beliefs, KL direction, and numerical convergence diagnostics. R-squared
#'   beliefs add coefficient/covariance references, numerical verification and
#'   observed-variance diagnostics that are never fitting constraints.
#' @param ... Reserved.
#' @return A summary.flexbayes_priors object, displayed by its print method.
#' @export
summary.flexbayes_priors = function(object, diagnostics = FALSE, ...) {
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) stop("diagnostics must be TRUE or FALSE.")
  structure(list(object = object, diagnostics = diagnostics), class = "summary.flexbayes_priors")
}

#' @export
print.summary.flexbayes_priors = function(x, ...) {
  if (x$diagnostics) prior_package_switch_summary(x$object)
  if (!is.null(x$object$joint)) return(prior_joint_summary(x))
  object = x$object
  diagnostics = x$diagnostics
  print(object)
  prior_r_squared_summary(object, diagnostics)
  for (name in names(object$beliefs)) {
    records = object$beliefs[[name]]
    if (!length(records) && is.null(object$fallbacks[[name]])) next
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
    outcome_records = Filter(function(r) prior_outcome_scale(r$scale), records)
    ratio_records = Filter(function(r) prior_ratio_scale(r$scale), outcome_records)
    if (length(ratio_records)) {
      prior_console_line(paste0(
        "Requested scale: ",
        ratio_records[[1]]$scale,
        "; multiplicative effect of a one-unit predictor increase, holding other predictors constant. This is not a probability ratio or an additive response change."
      ))
    }
    if (length(outcome_records) && !length(ratio_records)) {
      prior_console_line("Outcome-scale beliefs describe the expected outcome at the model-matrix centering reference, not individual observed outcomes.")
    }
    candidate = object$candidates[[name]]
    retained = !is.null(object$fallbacks[[name]])
    if (retained) prior_fallback_summary(object, name, diagnostics) else
      if (!is.null(candidate$message)) prior_console_line(paste0("Warning (", name, "): ", candidate$message))
    if (name == "sigma" && prior_has_sigma(object)) prior_scale_summary(object, diagnostics)
    if (diagnostics && length(ratio_records)) prior_ratio_diagnostics(object, name)
    if (diagnostics && identical(candidate$status, "contradictory")) {
      cat("  Candidate fitting skipped: mathematical contradiction.\n")
    } else if (diagnostics) {
      cat(
        "  ",
        prior_candidate_label(candidate$family),
        " candidate: ",
        if (retained) "not selected (fallback retained)" else if (isTRUE(candidate$valid)) "accepted" else "rejected",
        "\n",
        sep = ""
      )
      if (length(candidate$alternatives)) {
        for (alternative in candidate$alternatives) {
          prior_console_line(paste0(
            "  ",
            alternative$family,
            ": ",
            if (isTRUE(alternative$valid)) paste0("acceptable, mathematical rank ", alternative$rank) else paste0("rejected; ", alternative$message)
          ))
        }
        if (!is.null(candidate$mathematical)) prior_console_line(paste("Mathematical preference:", prior_format(candidate$mathematical)))
      }
      if (length(outcome_records)) {
        induced = if (outcome_records[[1]]$scale == "probability") "logit-normal / logit-Student's t" else "lognormal / exponentiated Student's t"
        prior_console_line(paste("Induced candidate families:", induced))
        if (!is.null(candidate$numerical)) prior_console_line(candidate$numerical)
        if (!is.null(candidate$convergence)) prior_console_line(paste("Optimizer convergence code:", candidate$convergence))
      }
      if (!is.null(candidate$score)) cat("  Score: ", format(candidate$score, digits = 4), "\n", sep = "")
      if (!is.null(candidate$achieved)) {
        details = data.frame(
          Belief = seq_along(records),
          Target = vapply(records, `[[`, numeric(1), "target"),
          Achieved = candidate$achieved,
          Tolerance = vapply(records, `[[`, numeric(1), "tolerance"),
          Residual = candidate$residuals
        )
        print(details, row.names = FALSE, digits = 4)
        prior_console_line("Residuals are divided by each belief's tolerance; acceptable absolute values are <= 1.")
      } else prior_console_line("Constraint diagnostics unavailable: candidate rejected before evaluation.")
    }
  }
  if (diagnostics) {
    cat("\n")
    prior_console_line(prior_default_description(object))
    if ("(Intercept)" %in% names(object$selected)) {
      prior_console_line("Centered intercept reference: model-matrix column means (including contrast and interaction columns).")
      reference = object$model$reference
      if (length(reference)) print(
        data.frame(Column = names(reference), Reference = unname(reference)),
        row.names = FALSE,
        digits = 4
      ) else cat("Intercept-only model; no predictor columns.\n")
    } else cat("No intercept; predictors are not centered.\n")
  }
  invisible(x)
}

#' Plot selected and candidate prior densities
#'
#' `plot(priors)` shows the selected distribution for each parameter. Supply
#' `parameter` to compare acceptable candidates for one coefficient, such as
#' `plot(priors, parameter = "x")`. These are parameter densities; to examine
#' their implications for outcomes and relationships, plot a `check_prior()` result.
#' After unsuccessful elicitation, plots show the retained prior with a fallback
#' annotation, rather than displaying an unselected candidate as the current prior.
#' @param x Prior specification.
#' @param parameter Optional parameter name. When supplied, compare acceptable candidates.
#' @param scale Plot in original, standardized, probability, response, odds_ratio
#'   or rate_ratio units. With parameter = 'sigma', r_squared shows the induced
#'   numerical cumulative distribution under the empirical predictor distribution;
#'   ordinary sigma density plots are unchanged. Ratio plots use a positive log10 x-axis.
#'   Outcome scales require the matching GLM centered intercept; ratios require
#'   an eligible numeric slope. Transformed densities
#'   include the transformation Jacobian.
#' @param ... Reserved.
#' @return A ggplot object.
#' @seealso [plot.flexbayes_prior_predictive()], [check_prior()]
#' @examples
#' \dontrun{
#' priors = identify_priors(mpg ~ wt, data = mtcars)
#' plot(priors)
#' plot(priors, parameter = "wt")
#' }
#' @export
plot.flexbayes_priors = function(x, parameter = NULL, scale = "original", ...) {
  if (length(scale) != 1L ||
      !is.character(scale) ||
      is.na(scale) ||
      !scale %in% c("original", "standardized", "probability", "response", "odds_ratio", "rate_ratio", "r_squared")) {
    stop("scale must be original, standardized, probability, response, odds_ratio, rate_ratio or r_squared.")
  }
  if (scale == "r_squared") {
    if (!identical(parameter, "sigma") || !prior_has_sigma(x)) stop("Select parameter = 'sigma' for a Gaussian R-squared plot.")
    return(prior_r_squared_plot(x))
  }
  factor = 1
  if (prior_outcome_scale(scale)) {
    if (is.null(parameter) || length(parameter) != 1L) stop("Select one eligible parameter for a transformed-scale plot.")
    prior_outcome_validate(list(scale = scale, type = "mean", target = if (scale == "probability") .5 else 1), x, parameter)
  }
  if (scale == "standardized") {
    if (is.null(parameter) || length(parameter) != 1L) stop("Select one parameter for a standardized plot.")
    factor = prior_standardization_factor(x, parameter)
  }
  names = if (is.null(parameter)) names(x$selected) else parameter
  if (!length(names) || any(!names %in% names(x$selected))) stop("Unknown parameter.")
  rows = list()
  for (name in names) {
    distributions = list(Selected = x$selected[[name]])
    if (!is.null(parameter) && is.null(x$fallbacks[[name]])) {
      alternatives = Filter(function(c) isTRUE(c$valid), x$candidates[[name]]$alternatives)
      distributions = lapply(alternatives, `[[`, "distribution")
      if (!length(distributions)) distributions = list(Selected = x$selected[[name]])
    }
    distributions = Filter(Negate(is.null), distributions)
    if (!length(distributions)) stop("No selected prior for ", name)
    for (label in names(distributions)) {
      p = distributions[[label]]
      grid = prior_distribution_operation(p, "quantile", seq(.001, .999, length.out = 501))
      if (name == "sigma" && prior_has_sigma(x) && any(!is.finite(grid))) {
        stop("Sigma prior density exceeds the numerical plotting range; no density points were silently removed.")
      }
      density = prior_distribution_operation(p, "density", grid)
      if (name == "sigma" && prior_has_sigma(x) && any(!is.finite(density))) {
        stop("Sigma prior density exceeds the numerical plotting range; no density points were silently removed.")
      }
      if (prior_outcome_scale(scale)) {
        outcome = if (scale == "probability") stats::plogis(grid) else exp(grid)
        jacobian = if (scale == "probability") outcome * (1 - outcome) else outcome
        density = density / jacobian
        keep = is.finite(outcome) & is.finite(density) & jacobian > 0
        if (any(!keep) && prior_ratio_scale(scale)) stop("Ratio density exceeds the numerical plotting range; no density points were silently removed.")
        if (any(!keep)) warning("Some transformed density points reached numerical boundaries and were omitted from the plot.", call. = FALSE)
        grid = outcome[keep]
        density = density[keep]
        if (!length(grid)) stop("Transformed density cannot be displayed reliably at this numerical scale.")
      }
      selected = identical(p, x$selected[[name]])
      rows[[length(rows) + 1L]] = data.frame(value = grid / factor, density = density * factor, parameter = name,
        Candidate = paste0(label, if (!is.null(x$fallbacks[[name]]) && selected) " (retained fallback)" else if (selected) " (selected)" else " (mathematical only)"))
    }
  }
  result = ggplot2::ggplot(do.call(rbind, rows), ggplot2::aes(x = value, y = density, colour = Candidate)) + ggplot2::geom_line() + ggplot2::facet_wrap(~parameter, scales = "free") + ggplot2::labs(
    x = if (scale == "standardized") "Standardized coefficient" else if (scale == "probability") "Probability at centering reference" else if (scale == "response") "Expected count at centering reference" else if (scale == "odds_ratio") "Odds ratio (one-unit predictor increase)" else if (scale == "rate_ratio") "Rate ratio (one-unit predictor increase)" else if (x$model$family$family == "binomial") "Coefficient (log-odds)" else if (x$model$family$family == "poisson") "Coefficient (log-rate)" else "Parameter value",
    y = "Prior density"
  )
  if (prior_ratio_scale(scale)) result = result + ggplot2::scale_x_log10()
  retained = intersect(names, prior_fallback_names(x))
  if (length(retained)) result = result + ggplot2::labs(subtitle = paste("Retained fallback:", paste(retained, collapse = ", "), "- current beliefs not verified"))
  result
}
