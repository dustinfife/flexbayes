# Keep a compact recipe: outcome simulations still use the original prediction cases.
prior_model_grid_policy = function() list(numeric = 5L, x = 51L, low_cardinality = 10L,
  discrete = 9L, maximum = 8192L, minimum = 2L, prediction_batch = 100000L)

prior_grid_numeric_values = function(values, count, discrete) {
  values = sort(unique(values))
  if (discrete) return(values[prior_display_indices(length(values), count)])
  if (length(values) == 1L) return(values)
  fraction = seq(0, 1, length.out = count)
  min(values) * (1 - fraction) + max(values) * fraction
}

prior_model_grid_specification = function(object, newdata) {
  variables = all.vars(stats::delete.response(object$model$terms))
  policy = prior_model_grid_policy()
  result = list(source = if (is.null(newdata)) "automatic" else "explicit", predictors = variables,
    policy = policy, weighting = if (is.null(newdata)) "Equal hypothetical grid configurations" else "Equal supplied cases within displayed coordinates")
  if (!is.null(newdata)) return(result)
  data = object$model$data
  if (!all(variables %in% names(data))) {
    result$error = "Automatic model grids require named predictor columns in the original data."
    return(result)
  }
  result$values = result$x_values = result$types = stats::setNames(vector("list", length(variables)), variables)
  for (name in variables) {
    value = data[[name]]
    if (!is.null(dim(value)) || !(is.numeric(value) || is.factor(value) || is.character(value) || is.logical(value))) {
      result$error = paste("Automatic model grids do not support predictor class for", name, "; supply supported explicit newdata.")
      return(result)
    }
    if (anyNA(value) || (is.numeric(value) && any(!is.finite(value)))) {
      result$error = paste("Nonfinite or missing original predictor values for", name)
      return(result)
    }
    if (is.numeric(value)) {
      discrete = all(value == floor(value)) || length(unique(value)) <= policy$low_cardinality
      result$types[[name]] = if (discrete) "discrete numeric" else "continuous numeric"
      result$values[[name]] = prior_grid_numeric_values(value, if (discrete) policy$discrete else policy$numeric, discrete)
      result$x_values[[name]] = prior_grid_numeric_values(value, policy$x, discrete)
    } else {
      result$types[[name]] = "categorical"
      categories = if (is.factor(value)) levels(value) else sort(unique(value))
      representative = if (is.factor(value)) factor(categories, levels = categories, ordered = is.ordered(value)) else categories
      result$values[[name]] = result$x_values[[name]] = representative
    }
  }
  result
}

prior_model_grid_expand = function(specification, axis) {
  if (!is.null(specification$error)) stop(specification$error)
  values = specification$values
  if (axis %in% names(values)) values[[axis]] = specification$x_values[[axis]]
  original = lengths(values)
  counts = original
  policy = specification$policy

  # Preserve the x resolution and every categorical level. Reduce numeric
  # nuisance resolution symmetrically, never truncate a Cartesian product.
  while (prod(as.double(counts)) > policy$maximum) {
    eligible = which(names(counts) != axis & unlist(specification$types) != "categorical" & counts > policy$minimum)
    if (!length(eligible)) stop("Balanced model grid exceeds ", policy$maximum,
      " configurations even at minimum numeric resolution. Supply a deliberate, smaller newdata grid; no combinations were silently dropped.")
    index = eligible[which.max(counts[eligible])]
    counts[index] = counts[index] - 1L
  }
  reduced = names(counts)[counts < original]
  for (name in reduced) {
    values[[name]] = prior_grid_numeric_values(values[[name]], counts[name], specification$types[[name]] != "continuous numeric")
  }
  if (length(reduced)) warning("Balanced model grid reduced numeric resolution (",
    paste(paste0(reduced, ": ", original[reduced], " -> ", lengths(values)[reduced]), collapse = ", "),
    "); equally weighted averages approximate the reduced grid. Inspect the plot's prior_model_grid attribute or supply explicit newdata.", call. = FALSE)
  data = if (length(values)) do.call(expand.grid, c(values, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))) else data.frame(row.names = 1L)
  list(data = data, values = values, original_counts = original, counts = lengths(values), reduced = reduced)
}

prior_model_coordinate_keys = function(grid, cases, columns) {
  codes = lapply(columns, function(name) {
    combined = c(as.character(grid[[name]]), as.character(cases[[name]]))
    match(combined, unique(combined))
  })
  key = do.call(paste, c(codes, list(sep = ":")))
  list(grid = head(key, nrow(grid)), cases = tail(key, nrow(cases)))
}

prior_model_marginal_grid = function(prep, information) {
  columns = unique(unlist(prep$display_roles, use.names = FALSE))
  columns = columns[!is.na(columns) & nzchar(columns)]
  keys = prior_model_coordinate_keys(prep$grid, prep$plot_data, columns)
  if (anyDuplicated(keys$grid)) stop("Model preparation produced duplicate displayed prediction coordinates.")
  group = match(keys$cases, keys$grid)
  if (anyNA(group)) stop("Prepared predictor cases do not match the displayed grid.")
  groups = split(seq_len(nrow(prep$data)), factor(group, levels = seq_len(nrow(prep$grid))))
  if (any(!lengths(groups))) stop("Model preparation created prediction configurations absent from the supplied cases.")
  information$omitted = setdiff(prep$variables$predictors, all.vars(prep$formula)[-1])
  attr(prep$grid, "prior_marginal") = list(cases = prep$data, groups = unname(groups))
  prep$marginal = information
  prep
}

prior_model_grid_subset = function(grid, rows) {
  result = grid[rows, , drop = FALSE]
  marginal = attr(grid, "prior_marginal")
  if (is.null(marginal)) return(result)
  marginal$groups = marginal$groups[rows]
  attr(result, "prior_marginal") = marginal
  result
}

prior_model_marginal_expected = function(x, marginal, ids) {
  groups = marginal$groups
  rows = unlist(groups, use.names = FALSE)
  membership = rep(seq_along(groups), lengths(groups))
  result = matrix(0, length(ids), length(groups))
  batch = max(1L, floor(prior_model_grid_policy()$prediction_batch / length(ids)))
  for (start in seq(1, length(rows), by = batch)) {
    index = seq.int(start, min(length(rows), start + batch - 1L))
    values = tryCatch(prior_model_expected(x, marginal$cases[rows[index], , drop = FALSE], ids),
      error = prior_model_configuration_error)
    sums = t(rowsum(t(values), membership[index], reorder = FALSE))
    slots = unique(membership[index])
    result[, slots] = result[, slots, drop = FALSE] + sums
  }
  sweep(result, 2, lengths(groups), "/")
}

prior_model_configuration_error = function(error) {
  stop("Cannot evaluate model-plot prediction configurations: ", conditionMessage(error),
    " Supply explicit newdata restricting predictions to valid configurations.", call. = FALSE)
}
