#' Evaluate Bayesian hypotheses
#'
#' Calculate the posterior probability that a hypothesis about
#' regression coefficients is true.
#'
#' @param object A brmsfit or stanreg model.
#' @param hypothesis An inequality involving model coefficients.
#' @param plot Logical. Plot the posterior distribution?
#'
#' @return An object of class bayes_hypothesis.
#'
#' @examples
#' \dontrun{
#' hypothesize(model, advertising > 1)
#' hypothesize(model, advertising > 1, plot = TRUE)
#' }
#'
#' @export
hypothesize = function(object, hypothesis, plot = FALSE) {
  
  expression = substitute(hypothesis)
  
  result = bayes_evaluate_hypothesis(
    object,
    expression
  )
  
  if (plot) {
    result$plot = bayes_plot_hypothesis(result)
    print(result$plot)
  }
  
  result
}


#' @export
print.bayes_hypothesis = function(x, ...) {
  
  cat("\nBayesian Hypothesis\n")
  cat("-------------------\n")
  
  cat("Hypothesis:", x$hypothesis, "\n")
  
  cat(
    "Posterior probability:",
    sprintf("%.3f", x$probability),
    "\n"
  )
  
  cat(
    "Percentage:",
    sprintf("%.1f%%", 100 * x$probability),
    "\n\n"
  )
  
  invisible(x)
}


#' @export
print.bayes_hypothesis = function(x, ...) {
  
  cat("\nBayesian Hypothesis\n")
  cat("-------------------\n")
  
  cat("Hypothesis:", x$hypothesis, "\n")
  
  cat(
    "Posterior probability:",
    sprintf("%.3f", x$probability),
    "\n"
  )
  
  cat(
    "Percentage:",
    sprintf("%.1f%%", 100 * x$probability),
    "\n\n"
  )
  
  invisible(x)
}