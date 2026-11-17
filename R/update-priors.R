#' Compose parameter-specific prior beliefs
#'
#' Construct a reusable specification and add it to a prior object with `+`.
#' Specifications do not fit distributions until applied. Addition delegates to
#' the same machinery as update(), preserving existing beliefs and diagnostics.
#' Chain additions as `priors + update_priors(...) + update_priors(...)`.
#' Addition requires a prior object on the left and a specification on the right;
#' specifications are not combined independently of a prior object.
#'
#' @param parameter One exact parameter name, as displayed by the prior object.
#' @param ... One or more objects constructed by belief().
#' @return update_priors() returns a flexbayes_prior_update specification.
#'   Adding it to a flexbayes_priors object returns a revised prior object.
#'   Neither input is modified. Unsupported or incompatible beliefs retain the
#'   same unresolved status as with update(); no constraints are discarded.
#' @examples
#' change = update_priors("x", belief(mode = 3), belief(sd = 1))
#' \dontrun{
#' priors = identify_priors(y ~ x + z, data = d) +
#'   change +
#'   update_priors("z", belief(mode = -1))
#' }
#' @export
update_priors = function(parameter, ...) {
  if (!is.character(parameter) || length(parameter) != 1L || is.na(parameter) || !nzchar(parameter)) stop("parameter must be one nonempty parameter name.")
  beliefs = list(...)
  if (!length(beliefs)) stop("Supply at least one belief() object.")
  if (any(!vapply(beliefs, inherits, logical(1), "prior_belief"))) stop("All additions must be belief objects.")
  structure(list(parameter = parameter, beliefs = beliefs), class = "flexbayes_prior_update")
}

#' @export
`+.flexbayes_priors` = function(e1, e2) {
  if (missing(e2)) stop("Use priors + update_priors(...); unary + is not supported.")
  if (!inherits(e1, "flexbayes_priors") || !inherits(e2, "flexbayes_prior_update")) stop("Use a flexbayes_priors object on the left and an update_priors() specification on the right.")
  # Unname the list so a belief argument cannot masquerade as parameter/remove.
  do.call(stats::update, c(list(object = e1, parameter = e2$parameter), unname(e2$beliefs)))
}

# Use the identical method for right-hand dispatch as well, so unsupported
# operands receive the same error and our two classes cannot conflict in Ops.
#' @export
`+.flexbayes_prior_update` = `+.flexbayes_priors`
