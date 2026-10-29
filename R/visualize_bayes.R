#' Visualize Bayesian Regression Models
#'
#' Visualizes Bayesian regression models using flexplot's existing
#' multivariate plotting system.
#'
#' @param object A brmsfit or stanreg model.
#' @param plot Which plot to display. Currently only "model" is supported.
#' @param formula An optional flexplot formula.
#' @param posterior_lines Number of sampled posterior regression lines (default 0).
#' @param seed Optional seed for selecting posterior draws.
#' @param ... Additional arguments passed to flexplot.
#' 
#' @exportS3Method flexplot::visualize brmsfit
#' @importFrom flexplot visualize
NULL
#'
#' @export
visualize.brmsfit = function(object, plot = "model", formula = NULL, ...) {
  bayes_visualize(object, plot, formula, ...)
}


#' @rdname visualize.brmsfit
#' @export
#' @exportS3Method flexplot::visualize stanreg
visualize.stanreg = function(object, plot = "model", formula = NULL, ...) {
  bayes_visualize(object, plot, formula, ...)
}


bayes_visualize = function(object, plot = "model", formula = NULL,
                            posterior_lines = 0, seed = NULL, ...) {
  if (plot != "model") {
    stop("Currently, Bayesian visualize() only supports plot = 'model'.")
  }

  if (length(posterior_lines) != 1L || !is.numeric(posterior_lines) ||
      is.na(posterior_lines) || !is.finite(posterior_lines) ||
      posterior_lines < 0 || posterior_lines != floor(posterior_lines)) {
    stop("posterior_lines must be a nonnegative integer.")
  }

  data = bayes_visualize_data(object)
  if (is.null(formula)) formula = bayes_visualize_formula(object, data)

  grid = flexplot:::generate_predictors(
    data = data, formula = formula, model = object, ...
  )
  draws = bayes_visualize_epred(object, grid)
  predictions = grid
  predictions$prediction = colMeans(draws)
  predictions$model = "Posterior mean"

  p = flexplot::flexplot(
    formula,
    data = data,
    prediction = predictions,
    suppress_smooth = TRUE,
    ...
  )

  if (posterior_lines == 0) return(p)

  lines = bayes_visualize_posterior_lines(
    draws = draws, grid = grid, formula = formula, data = data,
    x = rlang::as_name(p$mapping$x),
    posterior_lines = posterior_lines, seed = seed, ...
  )

  x = rlang::as_name(p$mapping$x)
  geometry = if (is.numeric(data[[x]])) ggplot2::geom_line else ggplot2::geom_point
  p + geometry(
    data = lines,
    mapping = ggplot2::aes(x = .data[[x]], y = .data[["prediction"]],
                           group = .data[[".posterior_group"]]),
    inherit.aes = FALSE, alpha = 0.15, colour = "gray50"
  )
}
