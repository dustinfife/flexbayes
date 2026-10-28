#' Bayesian regression estimates
#'
#' @param object A fitted brmsfit or stanreg model.
#' @param mc Retained for compatibility with flexplot.
#' @return Bayesian regression estimates.
#' @export
estimates.brmsfit = function(object, mc = TRUE) {
  
  model = bayes_extract_model(object)
  
  bayes_assemble_estimates(model)
  
}


#' @rdname estimates.brmsfit
#' @export
estimates.stanreg = function(object, mc = TRUE) {
  
  model = bayes_extract_model(object)
  
  bayes_assemble_estimates(model)
  
}
#' @export
print.bayes_estimates = function(x, ...) {
  
  cat("Bayesian regression estimates (95% credible intervals)\n\n")
  
  cat("Model R-squared:\n")
  print(x$r.squared, row.names = FALSE)
  
  cat("\nRegression coefficients:\n")
  print(x$coefficients, row.names = FALSE)
  
  if (!is.null(x$semi.p)) {
    cat("\nSemi-partial R-squared by term:\n")
    print(x$semi.p, row.names = FALSE)
  }
  
  if (length(x$factor.summary)) {
    for (v in names(x$factor.summary)) {
      cat("\nMarginal means for", v, ":\n")
      print(x$factor.summary[[v]], row.names = FALSE)
      
      cat("\nPairwise differences for", v, ":\n")
      print(x$difference.matrix[[v]], row.names = FALSE)
    }
  }
  
  invisible(x)
}