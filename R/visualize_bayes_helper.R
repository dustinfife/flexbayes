bayes_visualize_errors = function(object, options) {
  
  if (!inherits(object, c("brmsfit", "stanreg"))) stop("Model must be a brmsfit or stanreg object.")
  
  family = stats::family(object)$family
  if (family != "gaussian") stop("Currently, visualize() supports Gaussian Bayesian models only.")
  
  if (!is.null(options$posterior_lines)) {
    n = options$posterior_lines
    if (length(n) != 1 || !is.numeric(n) || is.na(n) || !is.finite(n) ||
        n < 0 || n != floor(n)) stop("posterior_lines must be a nonnegative integer.")
  }
  
  return(invisible(NULL))
}


bayes_visualize_data = function(object) {
  
  if (inherits(object, "brmsfit")) {
    return(as.data.frame(object$data))
  }
  
  if (inherits(object, "stanreg")) {
    return(as.data.frame(stats::model.frame(object)))
  }
  
  stop("Unsupported model class.")
}

bayes_model_formula = function(object) {
  
  model_formula = stats::formula(object)
  
  if (inherits(model_formula, "brmsformula")) {
    return(model_formula$formula)
  }
  
  model_formula
}

bayes_visualize_formula = function(object, data) {
  
  model_formula = bayes_model_formula(object)
  
  response = all.vars(model_formula)[1]
  
  predictors = flexplot:::get_predictors(object)
  predictors = flexplot:::remove_nonlinear_terms(predictors)
  predictors = unique(predictors)
  
  flexplot:::make_flexplot_formula(
    predictors = predictors,
    outcome = response,
    data = data
  )
}

bayes_visualize_epred = function(object, grid) {
  
  if (inherits(object, "brmsfit")) {
    return(brms::posterior_epred(
      object,
      newdata = grid,
      re_formula = NA
    ))
  }
  
  if (inherits(object, "stanreg")) {
    return(rstanarm::posterior_epred(
      object,
      newdata = grid,
      re.form = NA
    ))
  }
  
  stop("Unsupported model class.")
}

bayes_visualize_setup = function(object, data, options) {
  
  formula = options$formula
  if (is.null(formula)) formula = stats::formula(object)
  
  variables = all.vars(formula)
  response = all.vars(stats::formula(object))[1]
  predictors = setdiff(variables, response)
  
  if (length(predictors) == 0) stop("No predictors identified.")
  if (!all(predictors %in% names(data))) stop("Formula contains unknown predictors.")
  
  numeric_predictors = predictors[vapply(data[predictors], is.numeric, logical(1))]
  categorical_predictors = setdiff(predictors, numeric_predictors)
  
  if (length(numeric_predictors) > 1) {
    stop("This version supports one numeric predictor at a time. Supply a simpler formula.")
  }
  
  x = if (length(numeric_predictors)) numeric_predictors[1] else predictors[1]
  
  return(list(
    data = data,
    response = response,
    x = x,
    categorical = setdiff(categorical_predictors, x),
    formula = formula
  ))
}

bayes_visualize_grid = function(object, setup) {
  
  data = setup$data
  predictors = all.vars(stats::delete.response(stats::terms(object)))
  predictors = intersect(predictors, names(data))
  
  grid = data[rep(1, 100), predictors, drop = FALSE]
  
  for (name in predictors) {
    variable = data[[name]]
    
    if (is.numeric(variable)) grid[[name]] = rep(mean(variable, na.rm = TRUE), 100)
    if (is.factor(variable)) grid[[name]] = factor(levels(variable)[1], levels = levels(variable))
    if (is.character(variable)) grid[[name]] = rep(variable[which(!is.na(variable))[1]], 100)
  }
  
  x = setup$x
  if (is.numeric(data[[x]])) {
    grid[[x]] = seq(min(data[[x]], na.rm = TRUE),
                    max(data[[x]], na.rm = TRUE), length.out = 100)
  } else {
    values = if (is.factor(data[[x]])) levels(data[[x]]) else unique(data[[x]])
    grid = grid[seq_along(values), , drop = FALSE]
    grid[[x]] = values
    if (is.factor(data[[x]])) grid[[x]] = factor(values, levels = levels(data[[x]]))
  }
  
  return(grid)
}

bayes_visualize_predictions = function(object, data, formula, ...) {
  
  grid = flexplot:::generate_predictors(
    data = data,
    formula = formula,
    model = object,
    ...
  )
  
  draws = bayes_visualize_epred(object, grid)
  
  predictions = grid
  predictions$prediction = colMeans(draws)
  predictions$model = "Posterior mean"
  
  predictions
}

bayes_visualize_prediction_summary = function(draws, grid) {
  
  grid$prediction = colMeans(draws)
  grid$lower = apply(draws, 2, stats::quantile, probs = 0.025)
  grid$upper = apply(draws, 2, stats::quantile, probs = 0.975)
  
  return(grid)
}

bayes_visualize_sample_lines = function(draws, grid, n = 30) {
  
  if (n == 0) return(NULL)
  
  selected = sample(seq_len(nrow(draws)), min(n, nrow(draws)))
  lines = grid[rep(seq_len(nrow(grid)), times = length(selected)), , drop = FALSE]
  
  lines$prediction = as.vector(t(draws[selected, , drop = FALSE]))
  lines$draw = rep(seq_along(selected), each = nrow(grid))
  
  return(lines)
}

bayes_visualize_model = function(object, setup, options) {
  
  grid = bayes_visualize_grid(object, setup)
  draws = bayes_visualize_predictions(object, grid)
  predictions = bayes_visualize_prediction_summary(draws, grid)
  
  x = setup$x
  y = setup$response
  
  p = ggplot2::ggplot(setup$data, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Bayesian Regression", x = x, y = y)
  
  if (is.numeric(setup$data[[x]])) {
    p = p +
      ggplot2::geom_ribbon(
        data = predictions,
        ggplot2::aes(x = .data[[x]], ymin = .data$lower, ymax = .data$upper),
        inherit.aes = FALSE, alpha = 0.2
      ) +
      ggplot2::geom_line(
        data = predictions,
        ggplot2::aes(x = .data[[x]], y = .data$prediction),
        inherit.aes = FALSE, linewidth = 1
      )
  } else {
    p = p +
      ggplot2::geom_pointrange(
        data = predictions,
        ggplot2::aes(x = .data[[x]], y = .data$prediction,
                     ymin = .data$lower, ymax = .data$upper),
        inherit.aes = FALSE, linewidth = 0.8
      )
  }
  
  n = if (is.null(options$posterior_lines)) 0 else options$posterior_lines
  
  if (n > 0 && is.numeric(setup$data[[x]])) {
    lines = bayes_visualize_sample_lines(draws, grid, n)
    
    p = p + ggplot2::geom_line(
      data = lines,
      ggplot2::aes(x = .data[[x]], y = .data$prediction, group = .data$draw),
      inherit.aes = FALSE, alpha = 0.15, linewidth = 0.4
    )
  }
  
  return(p)
}

bayes_visualize_residuals = function(object, data) {
  
  predictions = bayes_visualize_predictions(object, data)
  fitted = colMeans(predictions)
  
  response = all.vars(stats::formula(object))[1]
  observed = data[[response]]
  residuals = observed - fitted
  
  residual_data = data.frame(fitted = fitted, residuals = residuals)
  residual_data$absolute = abs(residual_data$residuals)
  
  histogram = ggplot2::ggplot(residual_data, ggplot2::aes(x = .data$residuals)) +
    ggplot2::geom_histogram(bins = 30, fill = "lightgray", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Residual Distribution", x = "Residuals")
  
  dependence = ggplot2::ggplot(
    residual_data, ggplot2::aes(x = .data$fitted, y = .data$residuals)
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(se = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Residual Dependence")
  
  scale_location = ggplot2::ggplot(
    residual_data, ggplot2::aes(x = .data$fitted, y = .data$absolute)
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(se = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Scale-Location", y = "Absolute Residuals")
  
  return(cowplot::plot_grid(histogram, dependence, scale_location, ncol = 3))
}