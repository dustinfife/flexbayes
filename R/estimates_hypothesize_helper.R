# ============================================================
# BAYESIAN HYPOTHESIS HELPERS
# ============================================================


# EXTRACT POSTERIOR DRAWS -------------------------------------

bayes_hypothesis_draws = function(object) {
  
  coef_and_param_names = bayes_hypothesis_errors(object)
  parameter_names = coef_and_param_names$parameter_names
  coefficients    = coef_and_param_names$coefficients
  draws           = coef_and_param_names$draws

  result = draws[ , parameter_names, drop = FALSE]
  
  # Allow users to reference coefficients without
  # package-specific prefixes.
  colnames(result) = coefficients
  
  # Normalize the intercept name.
  which_are_intercept = colnames(result) == "(Intercept)"
  colnames(result)[which_are_intercept] = "Intercept"
  
  result
}

bayes_hypothesis_errors = function(object) {
  
  if (!inherits(object, c("brmsfit", "stanreg"))) stop("Model must be a brmsfit or stanreg object.")
  
  draws = as.matrix(object)
  
  if (inherits(object, "brmsfit")) {
    
    # brms prefixes population-level coefficients with b_.
    coefficients = rownames(brms::fixef(object))
    parameter_names = paste0("b_", coefficients)
    
  } else {
    
    # rstanarm uses the original coefficient names.
    coefficients = names(stats::coef(object))
    parameter_names = coefficients
    
  }
  
  if (!all(parameter_names %in% colnames(draws))) stop("Could not identify all regression coefficients.")
  
  return(list(
    draws = draws,
    coefficients = coefficients,
    parameter_names = parameter_names
  ))
}

# EVALUATE HYPOTHESIS -----------------------------------------

bayes_evaluate_hypothesis = function(object, expression) {
  
  bayes_hypothesis_expression_errors(expression)
  
  operator = as.character(expression[[1]])
  draws = bayes_hypothesis_draws(object)
  
  left = bayes_evaluate_expression(expression[[2]], draws)
  right = bayes_evaluate_expression(expression[[3]], draws)
  
  if (any(!is.finite(left)) || any(!is.finite(right))) stop("Hypothesis produced non-finite values.")
  
  satisfied = switch(
    operator,
    ">" = left > right,
    "<" = left < right,
    ">=" = left >= right,
    "<=" = left <= right
  )
  
  plotting = bayes_hypothesis_plot_parameters(left, right, expression)
  
  result = list(
    hypothesis = paste(deparse(expression), collapse = ""),
    probability = mean(satisfied),
    draws = left - right,
    satisfied = satisfied,
    operator = operator,
    plot_draws = plotting$plot_draws,
    threshold = plotting$threshold,
    x_label = plotting$x_label
  )
  
  class(result) = "bayes_hypothesis"
  
  return(result)
}

bayes_hypothesis_expression_errors = function(expression) {
  
  allowed = c(">", "<", ">=", "<=")
  error = "Hypothesis must be an inequality, such as advertising > 1."
  
  if (!is.call(expression)) stop(error)
  if (length(expression) != 3) stop(error)
  if (!as.character(expression[[1]]) %in% allowed) stop(error)
  
  return(invisible(NULL))
}

bayes_evaluate_expression = function(expr, draws) {
  
  n = nrow(draws)
  
  # Numeric constants
  if (is.numeric(expr) && length(expr) == 1) return(rep(expr, n))

  
  # Model coefficients
  if (is.symbol(expr)) {
    
    name = as.character(expr)
    if (!name %in% colnames(draws)) {
      stop("Unknown model coefficient: ", name)
    }
    return(draws[, name])
  }
  
  # Reject unsupported expressions
  if (!is.call(expr)) {
    stop("Unsupported hypothesis expression.")
  }
  
  op = as.character(expr[[1]])
  if (op == "(") return(bayes_evaluate_expression(expr[[2]], draws))
  if (!op %in% c("+", "-", "*", "/")) stop("Unsupported operation: ", op)
  
  # Unary arithmetic (+x or -x)
  if (length(expr) == 2) {
    
    x = bayes_evaluate_expression(expr[[2]], draws)
    
    if (op == "-") return(-x)
    if (op == "+") return(x)
  }
  
  if (length(expr) != 3) stop("Invalid arithmetic expression.")
  
  # Binary arithmetic
  a = bayes_evaluate_expression(expr[[2]], draws)
  b = bayes_evaluate_expression(expr[[3]], draws)
  
  result = switch(
    op,
    "+" = a + b,
    "-" = a - b,
    "*" = a * b,
    "/" = a / b
  )
  
  return(result)
}

bayes_hypothesis_plot_parameters = function(left, right, expression) {
  
  # If the right side is constant, use its original scale.
  right_is_constant = is.numeric(expression[[3]]) && length(expression[[3]]) == 1
  
  if (right_is_constant) {
    
    return(list(
      plot_draws = left,
      threshold = right[1],
      x_label = paste(deparse(expression[[2]]), collapse = "")
    ))
  }
  
  # Otherwise, plot the difference.
  result = list(
    plot_draws = left - right,
    threshold = 0,
    x_label = paste0(
      "(",
      paste(deparse(expression[[2]]), collapse = ""),
      ") - (",
      paste(deparse(expression[[3]]), collapse = ""),
      ")"
    )
  )
  
  if (any(!is.finite(left)) || any(!is.finite(right))) stop("Hypothesis produced non-finite values.")
  
  return(result)
}

bayes_hypothesis_plot_data = function(result) {
  
  draws = result$plot_draws
  threshold = result$threshold
  
  if (length(unique(draws)) < 2) {
    stop("Cannot plot a constant hypothesis expression.")
  }
  
  density = stats::density(draws)
  
  plot_data = data.frame(
    x = density$x,
    y = density$y
  )
  
  greater = result$operator %in% c(">", ">=")
  
  # Interpolate density at the threshold so the
  # shaded region ends exactly at the reference line.
  
  threshold_y = stats::approx(
    x = plot_data$x,
    y = plot_data$y,
    xout = threshold,
    rule = 2
  )$y
  
  threshold_row = data.frame(
    x = threshold,
    y = threshold_y
  )
  
  if (greater) {
    
    shaded_data = plot_data[
      plot_data$x > threshold,
      ,
      drop = FALSE
    ]
    
  } else {
    
    shaded_data = plot_data[
      plot_data$x < threshold,
      ,
      drop = FALSE
    ]
  }
  
  shaded_data = rbind(
    shaded_data,
    threshold_row
  )
  
  shaded_data = shaded_data[
    order(shaded_data$x),
    ,
    drop = FALSE
  ]
  
  list(
    density = plot_data,
    shaded = shaded_data
  )
}

bayes_plot_hypothesis = function(result) {
  
  plot_data = bayes_hypothesis_plot_data(result)
  
  p = ggplot2::ggplot(
    plot_data$density,
    ggplot2::aes(x = x, y = y)
  ) +
    
    ggplot2::geom_area(
      data = plot_data$shaded,
      fill = "steelblue",
      alpha = 0.4
    ) +
    
    ggplot2::geom_line(
      linewidth = 1
    ) +
    
    ggplot2::geom_vline(
      xintercept = result$threshold,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    
    ggplot2::labs(
      title = paste(
        "Hypothesis:",
        result$hypothesis
      ),
      subtitle = sprintf(
        "Posterior probability = %.1f%%",
        100 * result$probability
      ),
      x = result$x_label,
      y = "Density"
    ) +
    
    ggplot2::theme_bw(
      base_size = 14
    )
  
  return(p)
}


# VISUALIZE HYPOTHESIS ----------------------------------------

bayes_plot_hypothesis = function(result) {
  
  draws = result$draws
  
  if (length(unique(draws)) < 2) {
    stop("Cannot plot a constant hypothesis expression.")
  }
  
  # Calculate posterior density
  density = stats::density(draws)
  
  plot_data = data.frame(
    x = density$x,
    y = density$y
  )
  
  # Identify region satisfying the hypothesis
  greater = result$operator %in% c(">", ">=")
  
  plot_data$region = if (greater) {
    plot_data$x >= 0
  } else {
    plot_data$x <= 0
  }
  
  # Create shaded density region
  shaded_data = plot_data[plot_data$region, ]
  
  p = ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
    
    ggplot2::geom_area(
      data = shaded_data,
      fill = "steelblue",
      alpha = 0.4
    ) +
    
    ggplot2::geom_line(
      linewidth = 1
    ) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    
    ggplot2::labs(
      title = paste("Hypothesis:", result$hypothesis),
      subtitle = sprintf(
        "Posterior probability = %.1f%%",
        100 * result$probability
      ),
      x = result$x_label,
      y = "Density"
    ) +
    
    ggplot2::theme_bw(base_size = 14)
  
  return(p)
}