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

prior_plot_model = function(x, response, formula, draws, max_display, alpha, bins, breaks, labels, display) {
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
  if (nrow(prep$grid) * (if (display == "interval") 1 else length(ids)) > max_display) stop("Prediction grid exceeds max_display; use fewer draws or a simpler visualization formula.")
  draw_column = tail(make.unique(c(names(prep$grid), all.vars(prep$formula), ".prior_draw")), 1)
  data = prep$data
  if (response == "simulated") {
    data = data[rep(seq_len(nrow(data)), each = length(ids)), , drop = FALSE]
    data[[prep$variables$outcome]] = as.vector(x$simulated[ids, , drop = FALSE])
    data[[draw_column]] = rep(ids, times = nrow(prep$data))
    data$.prior_observation = rep(seq_len(nrow(prep$data)), each = length(ids))
  }
  if (display == "draws") {
    values = prior_model_expected(x, prep$grid, ids)
    prediction = prior_model_draw_frame(prep$grid, values, ids, draw_column)
    p = prior_model_render(prep, data, prediction, draw_column, alpha, response)
  } else {
    summary = prior_model_summary(x, prep$grid, max_display)
    summary[[draw_column]] = 1L
    template = prior_model_render(prep, data, summary, draw_column, 1, response)
    overlay = prior_model_summary_layers(template, prep)
    if (display == "interval") {
      p = template
      p$layers = p$layers[!vapply(p$layers, prior_model_summary_layer, logical(1))]
    } else {
      values = prior_model_expected(x, prep$grid, ids)
      prediction = prior_model_draw_frame(prep$grid, values, ids, draw_column)
      p = prior_model_render(prep, data, prediction, draw_column, alpha, response)
    }
    # Bands go behind curves, medians remain prominent in front of them.
    p$layers = c(overlay$interval, p$layers, overlay$median)
  }
  p = prior_model_observed_range(p, x)
  p + ggplot2::theme_bw() + ggplot2::labs(
    title = if (response == "expected") "Expected relationship based on priors" else "Simulated outcomes based on priors",
    y = prior_plot_outcome_name(x), subtitle = NULL,
    caption = switch(display,
      draws = prior_plot_draw_note(length(ids)),
      interval = paste("Median and central 90% prior interval;", nrow(x$coefficient_draws), "prior draws."),
      both = paste("Central 90% prior interval from", nrow(x$coefficient_draws), "draws;", length(ids), "individual draws shown.")))
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

prior_model_render = function(prep, data, prediction, draw_column, alpha, response) {
  flexplot::flexplot(prep$formula, data, prediction = prediction,
    prediction.draw = draw_column, prediction.alpha = alpha,
    bins = prep$binning$bins, breaks = prep$binning$breaks, labels = prep$binning$labels,
    raw.data = response == "simulated", sample = if (response == "expected") 0 else Inf,
    suppress_smooth = TRUE, se = FALSE, jitter = FALSE, alpha = alpha)
}

prior_model_summary = function(x, grid, max_display) {
  if (any(c(".prior_lower", ".prior_upper") %in% names(grid))) stop("Predictor names .prior_lower and .prior_upper are reserved for model interval metadata.")
  ids = seq_len(nrow(x$coefficient_draws))
  # Bound temporary prediction matrices without subsampling parameter draws.
  batch = max(1, floor(max_display / length(ids)))
  quantiles = matrix(NA_real_, nrow(grid), 3)
  for (start in seq(1, nrow(grid), by = batch)) {
    rows = seq.int(start, min(nrow(grid), start + batch - 1))
    values = prior_model_expected(x, grid[rows, , drop = FALSE], ids)
    quantiles[rows, ] = t(apply(values, 2, stats::quantile, probs = c(.05, .5, .95), names = FALSE))
  }
  grid$prediction = quantiles[, 2]
  grid$.prior_lower = quantiles[, 1]
  grid$.prior_upper = quantiles[, 3]
  grid
}

prior_model_summary_layer = function(layer) {
  is.data.frame(layer$data) && all(c(".prior_lower", ".prior_upper") %in% names(layer$data)) &&
    (inherits(layer$geom, "GeomLine") || inherits(layer$geom, "GeomPoint"))
}

prior_model_summary_layers = function(template, prep) {
  interval = median = list()
  numeric_x = is.numeric(prep$plot_data[[prep$display_roles$x]])
  # Reuse returned layer mappings, binned columns and position objects. This
  # deliberately does not inspect private path names or recreate dodge slots.
  # If flexplot changes its returned layer contract, fail rather than misplace
  # intervals. Rendering tests cover this small object-structure dependency.
  for (layer in Filter(prior_model_summary_layer, template$layers)) {
    if (!nrow(layer$data)) next
    line = inherits(layer$geom, "GeomLine")
    if (!numeric_x && line) next
    mapping = layer$mapping
    bounds = mapping
    bounds$ymin = ggplot2::aes(ymin = .prior_lower)$ymin
    bounds$ymax = ggplot2::aes(ymax = .prior_upper)$ymax
    bounds$shape = NULL
    if (numeric_x && line) {
      bounds$fill = bounds$colour
      bounds$colour = bounds$linetype = NULL
      interval[[length(interval) + 1L]] = prior_model_summary_geom(ggplot2::geom_ribbon,
        layer, bounds, aesthetic = "fill", alpha = .15)
      median[[length(median) + 1L]] = prior_model_summary_geom(ggplot2::geom_line,
        layer, mapping, legend = TRUE, linewidth = .9)
    } else {
      bounds$linetype = NULL
      interval[[length(interval) + 1L]] = prior_model_summary_geom(ggplot2::geom_linerange,
        layer, bounds, linewidth = .6)
      median[[length(median) + 1L]] = prior_model_summary_geom(ggplot2::geom_point,
        layer, mapping, legend = TRUE, size = 2)
    }
  }
  if (!length(median)) stop("Unsupported flexplot prediction-layer structure for prior intervals; use display = 'draws'.")
  list(interval = interval, median = median)
}

prior_model_summary_geom = function(geom, layer, mapping, aesthetic = "colour", legend = FALSE, ...) {
  args = list(data = layer$data, mapping = mapping, position = layer$position,
    inherit.aes = FALSE, show.legend = legend, ...)
  if (!is.null(layer$aes_params$colour)) args[[aesthetic]] = layer$aes_params$colour
  do.call(geom, args)
}
