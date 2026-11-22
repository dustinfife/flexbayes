#' Compose parameter-specific prior beliefs
#'
#' Construct a reusable specification and add it to a prior object with `+`.
#' Specifications do not fit distributions until applied. Addition delegates to
#' the same machinery as update(), preserving existing beliefs and diagnostics.
#' Expected-outcome specifications can be added with update_outcome(); these
#' reconcile jointly and may explicitly relax coefficient beliefs. Updating a
#' coefficient recomputes reconciliation from its coefficient-only baseline.
#' Chain updates as `priors + update_priors(...) + update_priors(...)`.
#' Use this interface to add, remove or replace beliefs. Removal positions refer
#' to the existing list displayed by summary(), not internal belief identifiers.
#' All removals happen simultaneously before additions; surviving beliefs keep
#' their order and metadata, and new beliefs are appended in the supplied order.
#' Bounds are checked when the specification is applied to a prior object.
#' Removing beliefs reuses candidate selection, reconciliation and fallback
#' recovery. Removing all beliefs restores the starting default when compatible
#' with the remaining model; existing backend and outcome constraints still apply.
#' Addition requires a prior object on the left and a specification on the right;
#' specifications are not combined independently of a prior object.
#'
#' @param parameter One exact parameter name, as displayed by the prior object.
#' @param ... Objects constructed by belief(). May be empty for removal-only updates.
#' @param remove Unique positive integer positions of existing beliefs to remove.
#'   Defaults to integer() (none). Positions must be in range when applied.
#' @return update_priors() returns a flexbayes_prior_update specification.
#'   Adding it to a flexbayes_priors object returns a revised prior object.
#'   Neither input is modified. Unsupported or incompatible beliefs retain the
#'   same unresolved status as with update(); no constraints are discarded.
#' @examples
#' change = update_priors("x", belief(mode = 3), belief(sd = 1))
#' removal = update_priors("x", remove = 2)
#' replacement = update_priors("x", belief(sd = .5), remove = 2)
#' \dontrun{
#' priors = identify_priors(y ~ x + z, data = d) +
#'   change +
#'   update_priors("z", belief(mode = -1))
#' }
#' @export
update_priors = function(parameter, ..., remove = integer()) {
  if (!is.character(parameter) || length(parameter) != 1L || is.na(parameter) || !nzchar(parameter)) stop("parameter must be one nonempty parameter name.")
  beliefs = list(...)
  prior_validate_removals(remove)
  if (!length(beliefs) && !length(remove)) stop("Supply at least one belief() object or removal index.")
  if (any(!vapply(beliefs, inherits, logical(1), "prior_belief"))) stop("All additions must be belief objects.")
  structure(list(parameter = parameter, beliefs = beliefs, remove = remove), class = "flexbayes_prior_update")
}

#' @export
`+.flexbayes_priors` = function(e1, e2) {
  if (missing(e2)) stop("Use priors + update_priors(...); unary + is not supported.")
  if (inherits(e1, "flexbayes_priors") && inherits(e2, "flexbayes_outcome_update")) return(prior_apply_outcome(e1, e2))
  if (!inherits(e1, "flexbayes_priors") || !inherits(e2, "flexbayes_prior_update")) {
    stop("Use a flexbayes_priors object on the left and an update_priors() specification on the right.")
  }

  # Unname the list so a belief argument cannot masquerade as parameter/remove.
  remove = if ("remove" %in% names(e2)) e2$remove else integer()
  do.call(stats::update, c(list(object = e1, parameter = e2$parameter, remove = remove), unname(e2$beliefs)))
}

# Use the identical method for right-hand dispatch as well, so unsupported
# operands receive the same error and our two classes cannot conflict in Ops.
#' @export
`+.flexbayes_prior_update` = `+.flexbayes_priors`

#' @export
`+.flexbayes_outcome_update` = `+.flexbayes_priors`

# A specification has no belief-list length yet; application checks its bounds.
prior_validate_removals = function(remove, count = NULL) {
  if (!is.numeric(remove) || is.complex(remove) || anyNA(remove) || any(!is.finite(remove)) ||
      any(remove != floor(remove)) || any(remove < 1) || anyDuplicated(remove)) {
    stop("Invalid removal numbers: use unique positive integer belief positions.")
  }
  if (!is.null(count) && any(remove > count)) stop(
    "Invalid removal numbers: positions exceed the parameter's existing belief list (", count, " beliefs).")
  invisible(remove)
}
