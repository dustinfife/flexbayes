# ============================================================
# BAYESIAN ESTIMATES: HELPER FUNCTIONS
#
# Internal functions supporting estimates() for Bayesian
# regression models fitted with brms or rstanarm.
#
# Currently supports Gaussian linear regression.
# ============================================================


# MODEL EXTRACTION --------------------------------------------

bayes_extract_model = function(object) {
  
  package = if (inherits(object, "brmsfit")) {
    "brms"
  } else if (inherits(object, "stanreg")) {
    "rstanarm"
  } else {
    stop("Unsupported Bayesian model class.")
  }
  
  family = stats::family(object)$family
  
  if (family != "gaussian") {
    stop(
      "Bayesian estimates currently supports Gaussian models only."
    )
  }
  
  dat = insight::get_data(object)
  
  f = stats::formula(object)
  
  if (inherits(f, "brmsformula")) {
    f = f$formula
  }
  
  tt = stats::terms(f)
  
  response = all.vars(f[[2]])
  
  if (length(response) != 1 ||
      !response %in% names(dat)) {
    stop("A single observed response variable is required.")
  }
  
  model_variables = all.vars(tt)
  
  if (anyNA(dat[, model_variables, drop = FALSE])) {
    stop("Missing values in model variables are not supported.")
  }
  
  # This version supports ordinary Gaussian regression.
  # Random-effects terms require additional handling.
  
  if (length(insight::find_random(object)) > 0) {
    stop(
      "Mixed models are not yet supported by estimates()."
    )
  }
  
  draws = as.matrix(object)
  
  if (package == "brms") {
    
    coefficient_names = rownames(brms::fixef(object))
    
    draw_names = paste0("b_", coefficient_names)
    
  } else {
    
    coefficient_names = names(stats::coef(object))
    
    draw_names = coefficient_names
    
  }
  
  if (!all(draw_names %in% colnames(draws))) {
    stop("Could not identify fixed-effect coefficient draws.")
  }
  
  X = stats::model.matrix(
    stats::delete.response(tt),
    data = dat
  )
  
  predictors = all.vars(stats::delete.response(tt))
  
  factors = intersect(
    predictors,
    names(dat)[vapply(dat, is.factor, logical(1))]
  )
  
  numbers = intersect(
    predictors,
    names(dat)[vapply(dat, is.numeric, logical(1))]
  )
  
  model = list(
    object = object,
    package = package,
    family = family,
    formula = f,
    terms = tt,
    data = dat,
    response = response,
    draws = draws,
    coefficient_names = coefficient_names,
    draw_names = draw_names,
    X = X,
    factors = factors,
    numbers = numbers
  )
  
  # Calculate posterior predictions once.
  model$predictions = bayes_epred(model)
  
  # Calculate factor-level predictions once.
  model$factor_draws = lapply(
    model$factors,
    function(variable) {
      bayes_factor_draws(model, variable)
    }
  )
  
  names(model$factor_draws) = model$factors
  
  model
  
}



# POSTERIOR SUMMARIES -----------------------------------------

bayes_summarize = function(x, name = "") {
  
  x = as.numeric(x)
  
  q = stats::quantile(
    x,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  data.frame(
    parameter = name,
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    lower = q[1],
    upper = q[2],
    row.names = NULL,
    check.names = FALSE
  )
  
}


# POSTERIOR PREDICTIONS ---------------------------------------

bayes_epred = function(model, newdata = NULL) {
  
  if (is.null(newdata)) {
    newdata = model$data
  }
  
  if (model$package == "brms") {
    
    predictions = brms::posterior_epred(
      model$object,
      newdata = newdata
    )
    
  } else {
    
    predictions = rstanarm::posterior_epred(
      model$object,
      newdata = newdata
    )
    
  }
  
  if (length(dim(predictions)) != 2) {
    stop("Expected a draws-by-observations prediction matrix.")
  }
  
  if (nrow(predictions) != nrow(model$draws)) {
    stop("Posterior predictions and coefficient draws do not align.")
  }
  
  predictions
  
}


# MODEL R-SQUARED ---------------------------------------------

bayes_r_squared = function(model) {
  
  predictions = model$predictions
  
  fitted_variance = apply(
    predictions,
    1,
    stats::var
  )
  
  sigma = model$draws[, "sigma"]
  
  residual_variance = sigma^2
  
  r2 = fitted_variance /
    (fitted_variance + residual_variance)
  
  bayes_summarize(r2, "R_squared")
  
}


# SEMI-PARTIAL R-SQUARED --------------------------------------

bayes_semi_partial = function(model) {
  
  tt = model$terms
  
  term_labels = attr(tt, "term.labels")
  
  if (length(term_labels) <= 1) {
    return(NULL)
  }
  
  X = model$X
  
  assignment = attr(X, "assign")
  
  predictions = model$predictions
  
  fitted_variance = apply(
    predictions,
    1,
    stats::var
  )
  
  sigma2 = model$draws[, "sigma"]^2
  
  total_variance = fitted_variance + sigma2
  
  results = lapply(
    seq_along(term_labels),
    function(j) {
      
      # Remove every design-matrix column belonging
      # to the current term.
      
      keep = assignment != j
      
      X_reduced = X[, keep, drop = FALSE]
      
      qr_reduced = qr(X_reduced)
      
      # For each posterior draw, project fitted values
      # onto the reduced model's column space.
      
      unique_variance = vapply(
        seq_len(nrow(predictions)),
        function(k) {
          
          residual_component = qr.resid(
            qr_reduced,
            predictions[k, ]
          )
          
          stats::var(residual_component)
          
        },
        numeric(1)
      )
      
      semi_r2 = unique_variance / total_variance
      
      bayes_summarize(
        semi_r2,
        term_labels[j]
      )
      
    }
  )
  
  result = do.call(rbind, results)
  
  rownames(result) = NULL
  
  result
  
}


# REGRESSION COEFFICIENTS -------------------------------------

bayes_coefficient_estimates = function(model) {
  
  results = lapply(
    seq_along(model$coefficient_names),
    function(i) {
      
      coefficient = model$draws[
        ,
        model$draw_names[i]
      ]
      
      bayes_summarize(
        coefficient,
        model$coefficient_names[i]
      )
      
    }
  )
  
  result = do.call(rbind, results)
  
  rownames(result) = NULL
  
  result
  
}


# NUMERIC REGRESSION ESTIMATES --------------------------------

bayes_regression_estimates = function(model) {
  
  coefficients = bayes_coefficient_estimates(model)
  
  # Preserve the existing convention of returning
  # regression coefficients without the intercept.
  #
  # Coefficients involving factors or interactions
  # retain their usual conditional interpretation.
  
  coefficients = coefficients[
    !coefficients$parameter %in% c(
      "Intercept",
      "(Intercept)"
    ),
    ,
    drop = FALSE
  ]
  
  rownames(coefficients) = NULL
  
  coefficients
  
}


# FACTOR MARGINAL MEANS ---------------------------------------

bayes_factor_draws = function(model, variable) {
  
  dat = model$data
  
  levels_v = levels(dat[[variable]])
  
  n_draws = nrow(model$draws)
  
  level_draws = matrix(
    NA_real_,
    nrow = n_draws,
    ncol = length(levels_v),
    dimnames = list(NULL, levels_v)
  )
  
  for (j in seq_along(levels_v)) {
    
    newdata = dat
    
    # Set every observation to the current factor level,
    # retaining the observed values of other predictors.
    
    newdata[[variable]] = factor(
      rep(levels_v[j], nrow(newdata)),
      levels = levels_v,
      ordered = is.ordered(dat[[variable]])
    )
    
    predictions = bayes_epred(
      model,
      newdata = newdata
    )
    
    # Average predictions across observations separately
    # for each posterior draw.
    
    level_draws[, j] = rowMeans(predictions)
    
  }
  
  level_draws
  
}


bayes_factor_estimates = function(model) {
  
  results = list()
  
  for (variable in model$factors) {
    
    draws = bayes_factor_draws(
      model,
      variable
    )
    
    summaries = lapply(
      seq_len(ncol(draws)),
      function(j) {
        
        bayes_summarize(
          draws[, j],
          colnames(draws)[j]
        )
        
      }
    )
    
    results[[variable]] = do.call(
      rbind,
      summaries
    )
    
    rownames(results[[variable]]) = NULL
    
  }
  
  results
  
}


# PAIRWISE FACTOR DIFFERENCES ---------------------------------

bayes_mean_differences = function(model) {
  
  results = list()
  
  for (variable in model$factors) {
    
    draws = bayes_factor_draws(
      model,
      variable
    )
    
    levels_v = colnames(draws)
    
    if (length(levels_v) < 2) {
      next
    }
    
    pairs = utils::combn(
      seq_along(levels_v),
      2
    )
    
    summaries = lapply(
      seq_len(ncol(pairs)),
      function(j) {
        
        a = pairs[1, j]
        b = pairs[2, j]
        
        # Difference calculated within each posterior draw,
        # preserving the posterior dependence between means.
        
        difference = draws[, b] - draws[, a]
        
        name = paste(
          levels_v[b],
          "-",
          levels_v[a]
        )
        
        bayes_summarize(
          difference,
          name
        )
        
      }
    )
    
    results[[variable]] = do.call(
      rbind,
      summaries
    )
    
    rownames(results[[variable]]) = NULL
    
  }
  
  results
  
}


# RESIDUAL STANDARD DEVIATION ---------------------------------

bayes_sigma = function(model) {
  
  if (!"sigma" %in% colnames(model$draws)) {
    return(NULL)
  }
  
  bayes_summarize(
    model$draws[, "sigma"],
    "sigma"
  )
  
}


# ASSEMBLE RESULTS --------------------------------------------

bayes_assemble_estimates = function(model) {
  
  result = list(
    r.squared = bayes_r_squared(model),
    semi.p = bayes_semi_partial(model),
    factor.summary = bayes_factor_estimates(model),
    difference.matrix = bayes_mean_differences(model),
    factors = model$factors,
    numbers.summary = bayes_regression_estimates(model),
    numbers = model$numbers,
    coefficients = bayes_coefficient_estimates(model),
    sigma = bayes_sigma(model),
    model.comparison = NULL
  )
  
  class(result) = "bayes_estimates"
  
  attr(result, "ci") = 0.95
  
  result
  
}