prior_display_integer = function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) stop(name, " must be a positive integer.")
  invisible(value)
}

prior_display_indices = function(total, count) unique(as.numeric(round(seq(1, total, length.out = min(total, count)))))

prior_display_values = function(values, max_display) {
  indices = prior_display_indices(length(values), max_display)
  data.frame(value = values[indices], .draw = (indices - 1) %% nrow(values) + 1,
    .observation = (indices - 1) %/% nrow(values) + 1)
}

prior_plot_outcome = function(x, response, max_display) {
  data = prior_display_values(x[[response]], max_display)
  if (any(!is.finite(data$value))) stop("Cannot plot nonfinite prior predictive values.")
  p = ggplot2::ggplot(data, ggplot2::aes(x = value))
  discrete = response == "simulated" && x$metadata$family$family %in% c("binomial", "poisson")
  if (discrete) p = p + ggplot2::geom_bar(ggplot2::aes(y = ggplot2::after_stat(prop), group = 1)) else p = p + ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30)
  outcome = prior_plot_outcome_name(x)
  p + ggplot2::theme_bw() + ggplot2::labs(
    title = paste(c("Hypothetical distribution of", if (response == "expected") "expected", outcome, "based on priors"), collapse = " "),
    x = if (response == "expected") paste("Expected", outcome, "(based on priors)") else paste("Simulated", outcome, "based on priors"),
    y = if (discrete) "Proportion" else "Density",
    caption = prior_plot_draw_note(nrow(x[[response]])))
}

prior_model_preparation = function(x, formula, bins, breaks, labels) {
  if (!"prepare_flexplot" %in% getNamespaceExports("flexplot") ||
      !all(c("prediction.draw", "prediction.alpha") %in% names(formals(flexplot::flexplot)))) stop("Install a flexplot version providing prepare_flexplot() and prediction.draw rendering.")
  metadata = x$metadata
  if (is.null(metadata$prediction_data) || is.null(metadata$terms)) stop("This check lacks model-plot metadata; rerun check_prior().")
  if (!is.symbol(metadata$formula[[2]])) stop("Model plots currently require an untransformed response variable.")
  data = metadata$prediction_data
  original_design = prior_model_grid_matrix(x, data)
  if (!isTRUE(all.equal(original_design, x$X))) stop("Model terms or predictor values changed since check_prior(); rerun the check.")
  outcome = as.character(metadata$formula[[2]])
  # New predictor-only cases have no observations to hide or overwrite. Supply
  # stored expected medians solely for the public preparer's outcome requirement.
  if (!outcome %in% names(data)) data[[outcome]] = apply(x$expected, 2, stats::median)
  prep = flexplot::prepare_flexplot(data, model = metadata$formula, formula = formula,
    bins = bins, breaks = breaks, labels = labels)
  if (!identical(prep$version, 1L)) stop("Unsupported flexplot preparation contract version.")
  if (!length(prep$variables$predictors)) stop("Model plots require a predictor; use an outcome or interval plot for intercept-only models.")
  if (!identical(as.integer(prep$rows), seq_len(nrow(data)))) stop("Model plot preparation omitted cases with missing values; supply complete prediction cases.")
  prep
}

prior_model_grid_matrix = function(x, grid) {
  terms = stats::delete.response(x$metadata$terms)
  # model.frame rebuilds factor levels; matrix contrasts are restored explicitly below.
  for (name in names(grid)) if (is.factor(grid[[name]])) attr(grid[[name]], "contrasts") = NULL
  frame = stats::model.frame(terms, grid, na.action = stats::na.fail, xlev = x$metadata$factor_levels)
  X = stats::model.matrix(terms, frame, contrasts.arg = x$metadata$contrasts)
  if (!identical(colnames(X), colnames(x$coefficient_draws)) || any(!is.finite(X))) stop("Visualization grid does not match the original coefficient parameterization.")
  X
}

prior_model_expected = function(x, grid, ids) {
  X = prior_model_grid_matrix(x, grid)
  linear = x$coefficient_draws[ids, , drop = FALSE] %*% t(X)
  if (any(!is.finite(linear))) stop("Nonfinite model-plot linear predictors.")
  family = x$metadata$family$family
  expected = switch(family, gaussian = linear, binomial = stats::plogis(linear), poisson = exp(linear))
  if (any(!is.finite(expected)) || (family == "poisson" && any(expected > 1e8))) stop("Model-plot predictions exceed numerical limits (Poisson rates must not exceed 1e8).")
  expected
}

prior_model_draw_frame = function(grid, values, ids, draw_column) {
  data = grid[rep(seq_len(nrow(grid)), each = length(ids)), , drop = FALSE]
  data[[draw_column]] = rep(ids, times = nrow(grid))
  data$prediction = as.vector(values)
  data
}

prior_plot_model = function(x, response, formula, draws, max_display, alpha, bins, breaks, labels) {
  rng = prior_plot_rng_state()
  on.exit(prior_restore_plot_rng(rng), add = TRUE)
  prior_display_integer(draws, "draws")
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha < 0 || alpha > 1) stop("prediction.alpha must be between zero and one.")
  prep = prior_model_preparation(x, formula, bins, breaks, labels)
  count = min(draws, nrow(x$coefficient_draws))
  if (response == "simulated") {
    if (nrow(prep$data) > max_display) stop("Too many prediction cases for simulated model display; increase max_display or use fewer newdata cases.")
    count = min(count, floor(max_display / nrow(prep$data)))
  }
  ids = prior_display_indices(nrow(x$coefficient_draws), count)
  if (nrow(prep$grid) * length(ids) > max_display) stop("Prediction grid exceeds max_display; use fewer draws or a simpler visualization formula.")
  values = prior_model_expected(x, prep$grid, ids)
  draw_column = tail(make.unique(c(names(prep$grid), all.vars(prep$formula), ".prior_draw")), 1)
  prediction = prior_model_draw_frame(prep$grid, values, ids, draw_column)
  data = prep$data
  if (response == "simulated") {
    data = data[rep(seq_len(nrow(data)), each = length(ids)), , drop = FALSE]
    data[[prep$variables$outcome]] = as.vector(x$simulated[ids, , drop = FALSE])
    data[[draw_column]] = rep(ids, times = nrow(prep$data))
    data$.prior_observation = rep(seq_len(nrow(prep$data)), each = length(ids))
  }
  p = flexplot::flexplot(prep$formula, data, prediction = prediction,
    prediction.draw = draw_column, prediction.alpha = alpha,
    bins = prep$binning$bins, breaks = prep$binning$breaks, labels = prep$binning$labels,
    raw.data = response == "simulated", sample = if (response == "expected") 0 else Inf, suppress_smooth = TRUE, se = FALSE,
    jitter = FALSE, alpha = alpha)
  p = prior_model_observed_range(p, x)
  p + ggplot2::theme_bw() + ggplot2::labs(
    title = if (response == "expected") "Expected relationship based on priors" else "Simulated outcomes based on priors",
    y = prior_plot_outcome_name(x), subtitle = NULL,
    caption = prior_plot_draw_note(length(ids)))
}

prior_plot_rng_state = function() {
  present = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(present = present, seed = if (present) get(".Random.seed", envir = .GlobalEnv))
}
prior_restore_plot_rng = function(state) {
  if (state$present) assign(".Random.seed", state$seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
}

prior_plot_outcome_name = function(x) {
  outcome = x$metadata$formula[[2]]
  if (is.symbol(outcome)) return(as.character(outcome))
  paste(deparse(outcome, width.cutoff = 500L), collapse = " ")
}

prior_plot_draw_note = function(count) paste("Based on", count, "prior draws.")

prior_model_observed_range = function(p, x) {
  observed = x$metadata$observed_response
  # Older checks retain training outcomes only when no newdata was supplied.
  if (is.null(observed) && identical(x$metadata$newdata, FALSE)) {
    observed = x$metadata$prediction_data[[prior_plot_outcome_name(x)]]
  }
  if (is.null(observed)) {
    warning("Original observed outcomes unavailable; rerun check_prior() to display the observed range.", call. = FALSE)
    return(p)
  }
  observed = observed[is.finite(observed)]
  if (!length(observed)) {
    warning("No finite original observed outcomes; observed range references omitted.", call. = FALSE)
    return(p)
  }
  reference = data.frame(.observed_value = range(observed),
    .observed_label = c("Observed minimum", "Observed maximum"), .observed_vjust = c(1.2, -.4))
  p + ggplot2::geom_hline(data = reference, ggplot2::aes(yintercept = .observed_value),
    linetype = "dashed", colour = "grey45", linewidth = .4, alpha = .65, show.legend = FALSE) +
    ggplot2::geom_text(data = reference,
      ggplot2::aes(x = Inf, y = .observed_value, label = .observed_label, vjust = .observed_vjust),
      inherit.aes = FALSE, hjust = 1.05, colour = "grey45", size = 3, show.legend = FALSE)
}
