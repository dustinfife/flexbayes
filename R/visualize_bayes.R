#' Visualize Bayesian Regression Models
#'
#' Visualizes Bayesian regression models using flexplot's existing
#' multivariate plotting system.
#'
#' @param object A brmsfit or stanreg model.
#' @param plot Which plot to display. Currently only "model" is supported.
#' @param formula An optional flexplot formula.
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


bayes_visualize = function(object, plot = "model", formula = NULL, ...) {
  if (plot != "model") {
    stop("Currently, Bayesian visualize() only supports plot = 'model'.")
  }

  data = bayes_visualize_data(object)
  if (is.null(formula)) formula = bayes_visualize_formula(object, data)

  predictions = bayes_visualize_predictions(object, data, formula, ...)

  flexplot::flexplot(
    formula,
    data = data,
    prediction = predictions,
    suppress_smooth = TRUE,
    ...
  )
}
