# Difference specifications reuse outcome validation and the same joint solver.
prior_difference_update = function(difference, records) {
  if (!is.list(difference) || length(difference) != 2L ||
      is.null(names(difference)) || anyDuplicated(names(difference)) ||
      !setequal(names(difference), c("from", "to"))) {
    stop("difference must be a named list containing exactly from and to configurations.")
  }
  prior_outcome_update(difference$from, records, "difference$from")
  prior_outcome_update(difference$to, records, "difference$to")
  list(difference = difference[c("from", "to")], beliefs = unname(records))
}

# Both endpoints pass through the original prediction matrix and centering.
prior_joint_configuration = function(object, at, variables, target) {
  if (!setequal(names(at), variables)) stop(target, " must specify exactly these predictors: ", paste(variables, collapse = ", "))
  data = if (length(variables)) as.data.frame(at) else data.frame(row.names = 1L)
  X = prior_prediction_matrix(object, data)
  slopes = setdiff(colnames(X), "(Intercept)")
  X[, slopes] = sweep(X[, slopes, drop = FALSE], 2, object$model$reference[slopes], "-")
  X
}

# A zero contrast is a point mass, so strict and inclusive comparisons differ.
prior_joint_constant_value = function(record) {
  if (record$type != "probability") return(0)
  as.numeric(switch(record$operator, `<` = 0 < record$threshold,
    `<=` = 0 <= record$threshold, `>` = 0 > record$threshold,
    `>=` = 0 >= record$threshold))
}

prior_joint_configuration_text = function(at) {
  if (!length(at)) return("intercept only")
  paste(names(at), vapply(at, prior_joint_summary_number, character(1)), sep = " = ", collapse = ", ")
}

prior_joint_summary_configuration = function(specification) {
  if (is.null(specification$difference)) {
    prior_console_line(paste("At:", prior_joint_configuration_text(specification$at)))
    return(invisible(NULL))
  }
  prior_console_line(paste("From:", prior_joint_configuration_text(specification$difference$from)))
  prior_console_line(paste("To:  ", prior_joint_configuration_text(specification$difference$to)))
  invisible(NULL)
}

# Keep the solver's constraint order independent of specification order, while
# retaining original configurations, IDs and belief positions for reporting.
prior_joint_constraint_key = function(record, problem) {
  paste(paste(sprintf("%.17g", problem$X[record$configuration, ]), collapse = "|"),
    prior_outcome_constraint_key(record), sep = "|")
}

prior_joint_ordered_problem = function(problem) {
  if (!isTRUE(problem$has_difference)) return(problem)
  keys = vapply(problem$outcome, prior_joint_constraint_key, character(1), problem = problem)
  indices = order(keys, method = "radix")
  problem$outcome = problem$outcome[indices]
  problem$units = problem$units[indices]
  problem
}
