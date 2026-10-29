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


# Apply flexplot's own prediction processing to each sampled draw.
# This preserves its binning, grouping and faceting conventions rather than
# maintaining a second implementation of those rules.
bayes_visualize_posterior_lines = function(draws, grid, formula, data, x,
                                           posterior_lines = 50, seed = NULL, ...) {
  if (!is.null(seed)) {
    selected = withr::with_seed(
      seed,
      sample.int(nrow(draws), min(posterior_lines, nrow(draws)))
    )
  } else {
    selected = sample.int(nrow(draws), min(posterior_lines, nrow(draws)))
  }
  
  processed = lapply(selected, function(i) {
    prediction = grid
    prediction$prediction = draws[i, ]
    prediction$model = "Posterior mean"
    
    plot = flexplot::flexplot(
      formula,
      data = data,
      prediction = prediction,
      suppress_smooth = TRUE,
      ...
    )
    
    # Identify layers containing processed prediction data.
    candidates = which(vapply(plot$layers, function(layer) {
      is.data.frame(layer$data) && "prediction" %in% names(layer$data)
    }, logical(1)))
    
    if (length(candidates) == 0L) {
      stop("Could not identify flexplot's prediction layer.")
    }
    
    # Categorical x-axes may produce both point and line layers.
    # Select the geometry appropriate for the predictor type.
    if (length(candidates) > 1L) {
      geom_class = if (is.numeric(data[[x]])) "GeomLine" else "GeomPoint"
      
      candidates = candidates[vapply(
        plot$layers[candidates],
        function(layer) inherits(layer$geom, geom_class),
        logical(1)
      )]
    }
    
    if (length(candidates) != 1L) {
      stop("Could not uniquely identify flexplot's prediction layer.")
    }
    
    result = plot$layers[[candidates]]$data
    result$.posterior_draw = i
    result
  })
  
  lines = do.call(rbind, processed)
  
  grouping = setdiff(
    names(lines),
    c(x, "prediction", "model", ".posterior_draw", ".posterior_group")
  )
  
  if (length(grouping)) {
    lines$.posterior_group = interaction(
      lines[c(".posterior_draw", grouping)],
      drop = TRUE,
      lex.order = TRUE
    )
  } else {
    lines$.posterior_group = factor(lines$.posterior_draw)
  }
  
  lines
}
