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
