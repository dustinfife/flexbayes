# Factors convert standardized slopes into original coefficient units.
prior_standardization_metadata = function(model) {
  if (model$family$family != "gaussian") return(list())
  result = list()
  terms = model$terms
  factors = attr(terms, "factors")
  symbols = Filter(is.symbol, as.list(attr(terms, "variables"))[-1])
  simple = vapply(symbols, as.character, character(1))
  for (variable in names(model$frame)) {
    if (!variable %in% names(model$data) || !variable %in% simple) next
    column = model$frame[[variable]]
    if (!is.numeric(column) || !is.null(dim(column))) next
    label = paste(deparse(as.name(variable), backtick = TRUE), collapse = "")
    term = match(label, attr(terms, "term.labels"))
    if (is.na(term)) next
    # A main effect involved in an interaction is outside this initial scope.
    row = match(label, rownames(factors))
    if (!is.na(row) && any(factors[row, ] > 0 & attr(terms, "order") > 1)) next
    indices = which(attr(model$X, "assign") == term)
    if (length(indices) != 1L) next
    sx = stats::sd(model$X[, indices])
    sy = stats::sd(stats::model.response(model$frame))
    result[[colnames(model$X)[indices]]] = list(predictor = variable,
      predictor_sd = sx, outcome_sd = sy, factor = sy / sx)
  }
  result
}

prior_standardization_factor = function(object, parameter) {
  if (!prior_has_sigma(object)) stop("Gaussian-style full standardization is not defined for these GLMs; specify beliefs on the linear predictor scale.")
  info = object$model$standardization[[parameter]]
  if (is.null(info)) stop("Standardized beliefs require an untransformed numeric slope without interactions; unsupported parameter: ", parameter)
  if (!is.finite(info$predictor_sd) || info$predictor_sd <= 0 ||
      !is.finite(info$outcome_sd) || info$outcome_sd <= 0 ||
      !is.finite(info$factor) || info$factor <= 0) stop("Standardization requires positive finite predictor and outcome SDs.")
  prior_validate_standardization_data(object)
  info$factor
}

prior_validate_standardization_data = function(object) {
  frame = stats::model.frame(object$model$formula, object$model$data, na.action = stats::na.fail)
  X = stats::model.matrix(attr(frame, "terms"), frame, contrasts.arg = object$model$contrasts)
  if (!isTRUE(all.equal(X, object$model$X)) ||
      !isTRUE(all.equal(stats::model.response(frame), stats::model.response(object$model$frame)))) stop("Model data changed since standardization; identify priors again.")
  invisible(object)
}

prior_validate_standardized_beliefs = function(object) {
  records = unlist(object$beliefs, recursive = FALSE)
  if (any(vapply(records, function(r) identical(r$scale, "standardized"), logical(1)))) prior_validate_standardization_data(object)
  invisible(object)
}

prior_convert_belief = function(record, object, parameter) {
  if (prior_outcome_scale(record$scale)) return(prior_outcome_validate(record, object, parameter))
  if (!identical(record$scale, "standardized")) return(record)
  factor = prior_standardization_factor(object, parameter)
  if (!is.null(record$entered)) {
    record$target = record$entered$target
    record$threshold = record$entered$threshold
    record$tolerance = record$entered$tolerance
  }
  record$entered = list(target = record$target, threshold = record$threshold,
    tolerance = record$tolerance, scale = record$scale)
  record$standardization_factor = factor
  if (record$type == "probability") record$threshold = record$threshold * factor else {
    record$target = record$target * factor
    record$tolerance = record$tolerance * factor
  }
  if (!is.finite(record$target) || (!is.null(record$threshold) && !is.finite(record$threshold)) || !is.finite(record$tolerance) || record$tolerance <= 0) stop("Standardized constraint exceeds numerical limits.")
  record
}
