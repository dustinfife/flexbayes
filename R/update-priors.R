#' Incorporate parameter, expected-outcome, or expected-difference beliefs
#'
#' Construct a reusable specification and add it to a prior object with `+`.
#' Specifications do not fit distributions until applied. Addition delegates to
#' the existing parameter-fitting or joint-reconciliation machinery, preserving
#' existing beliefs and diagnostics.
#' Supply exactly one target: parameter for an individual prior, outcome for
#' the conditional expected response, or difference for the change between
#' two expected responses. Outcome and difference
#' beliefs reconcile jointly and may explicitly relax coefficient beliefs. Updating a
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
#' If rstanarm restrictions block a belief update, flexbayes tries brms once and
#' switches only if the entire prior specification can be verified and exported.
#' This requires brms to be installed. Contradictory or numerically unresolved
#' beliefs retain the existing fallbacks and safeguards. There is no automatic
#' switch back; use switch_package() to request a package explicitly.
#'
#' @param parameter One exact parameter name, as displayed by the prior object.
#'   May also be supplied as the first unnamed argument.
#' @param outcome Named list with one finite numeric value for every predictor.
#'   Beliefs concern the conditional expected response, not individual noisy
#'   observations. For an intercept-only model use list(). Supply either this
#'   argument, parameter, or difference. Control argument names must match exactly.
#' @param difference Named list containing from and to, each a predictor
#'   configuration with the same requirements as outcome. Beliefs describe
#'   the conditional expected response at to minus that at from. Both responses
#'   share the same coefficients; residual variability is excluded. Exactly one
#'   of parameter, outcome, or difference must be supplied.
#' @param ... Objects constructed by belief(). May be empty for removal-only updates.
#' @param remove Unique positive integer positions of existing beliefs to remove.
#'   Defaults to integer() (none). Positions must be in range when applied.
#'   Only parameter targets support removal; do not supply remove with outcome
#'   or difference.
#' @section Expected-outcome reconciliation:
#' Initially supports Gaussian identity models with an intercept,
#' numeric predictors, and selected independent Normal coefficient priors.
#' Student's t coefficient priors, categorical predictors and offsets are
#' unsupported. Sigma priors can be elicited separately with update_priors();
#' expected-outcome reconciliation never adjusts sigma. The existing centered
#' intercept is retained; independence refers to that parameterization.
#'
#' Differences use the contrast between the two centered model-matrix rows.
#' The intercept and unchanged model-matrix columns cancel. Uncertainty uses
#' the shared coefficient vector, not the sum of independent endpoint variances.
#' Identical rows imply exactly zero, including strict versus inclusive
#' probability comparisons. Incompatible zero-difference beliefs are retained
#' and reported as mathematical incompatibility.
#'
#' Reconciliation has three priorities: satisfy outcome and difference beliefs within their
#' tolerances; minimize the largest coefficient-belief violation beyond its
#' tolerance in dimensionless coordinates; then minimize
#' KL(adjusted prior || baseline prior). Locations and SDs use baseline-SD
#' units; probabilities use probit (standard-Normal quantile) distances outside
#' their original probability tolerance bands. Numerical tolerances determine
#' satisfaction, not the relative importance of beliefs. Defaults enter
#' only the last tier. Compatibility checks first derive Normal moment
#' acceptance bounds and test nonnegative coefficient-variance restrictions.
#' A verified separating certificate establishes incompatibility; a failed
#' numerical search alone does not. Passing these necessary checks does not
#' prove full joint feasibility. All answers are retained when either check
#' or reconciliation fails. This is a deterministic constrained numerical search,
#' not a guarantee of global optimality. An unverified solution is never exported.
#' Original coefficient-only selections define the baseline on every update,
#' preventing successive outcome additions from moving the baseline.
#'
#' The default belief tolerance remains 1e-6, including for currency-valued
#' locations. For a price estimate stated to the nearest 10000, an explicit
#' tolerance=5000 may describe reporting precision. This is different from the
#' uncertainty expressed by a probability belief. Record such precision as part
#' of the interview, not as a way to force incompatible answers to fit.
#' Probability targets of 0 or 1 use their declared positive acceptance bands;
#' a finite Normal distribution cannot attain those probabilities exactly.
#'
#' The fixed-variance Mahalanobis projection initializes the numerical search.
#' Positive affine changes of simple numeric predictor units preserve the
#' centered independent-Normal structure. Equivalent targets, thresholds and
#' absolute tolerances give equivalent reconciliations up to numerical accuracy.
#' This does not extend to transformations mixing coefficient columns and
#' inducing covariance. Outcome unit changes must transform tolerances too.
#' General reconciliation may change both means and SDs. Underidentified problems
#' are resolved by the stated preference criterion, not by claiming identification.
#' Relaxed coefficient beliefs produce a warning and remain visible in
#' summary(..., diagnostics=TRUE). Incompatible or numerically unresolved outcome
#' specifications are retained but block check_prior(), as_rstanarm(), as_brms()
#' and fit() unless retained fallback priors are explicitly approved.
#' @return update_priors() returns a flexbayes_prior_update specification.
#'   Adding it to a flexbayes_priors object returns a revised prior object.
#'   Neither input is modified. Unsupported or incompatible beliefs retain the
#'   existing failure diagnostics and fallback safeguards; no constraints are discarded.
#' @examples
#' change = update_priors(parameter = "x", belief(mode = 3), belief(sd = 1))
#' removal = update_priors("x", remove = 2)
#' replacement = update_priors("x", belief(sd = .5), remove = 2)
#' expected = update_priors(outcome = list(x = 2), belief(mode = 10),
#'   belief(value < 8, probability = .05))
#' # An agent expects a 50000 price difference over 300 square feet and
#' # gives a 10% chance that the expected difference exceeds 90000.
#' # There is no need to convert these answers into a slope by hand.
#' houses = identify_priors(sale_price ~ age + sqft + bedrooms,
#'   data = homes_2024, package = "rstanarm") +
#'   update_priors(
#'     difference = list(
#'       from = list(sqft = 1700, age = 44, bedrooms = 3),
#'       to = list(sqft = 2000, age = 44, bedrooms = 3)),
#'     belief(mode = 50000),
#'     belief(value > 90000, probability = .10))
#' summary(houses)
#' \dontrun{
#' priors = identify_priors(y ~ x + z, data = d) +
#'   change +
#'   update_priors("z", belief(mode = -1))
#' }
#' @export
update_priors = function(..., parameter = NULL, outcome = NULL, difference = NULL, remove = integer()) {
  beliefs = list(...)
  has_parameter = !missing(parameter)
  has_outcome = !missing(outcome)
  has_difference = !missing(difference)

  # All control arguments match exactly after dots. Recognize only the leading
  # unnamed parameter shorthand; outcome calls keep every belief in dots.
  unnamed_first = length(beliefs) && (is.null(names(beliefs)) || !nzchar(names(beliefs)[1]))
  positional = unnamed_first && (is.character(beliefs[[1]]) ||
    (!has_parameter && !has_outcome && !has_difference && !inherits(beliefs[[1]], "prior_belief")))
  if (positional) {
    if (has_parameter) stop("Supply parameter only once, either positionally or by name.")
    parameter = beliefs[[1]]
    beliefs = beliefs[-1]
    has_parameter = TRUE
  }
  if (sum(has_parameter, has_outcome, has_difference) != 1L) stop("Supply exactly one of parameter, outcome, or difference, using exact argument names.")
  if (has_difference) {
    if (!missing(remove)) stop("remove applies only to parameter beliefs; difference removal is not supported.")
    specification = prior_difference_update(difference, beliefs)
    return(structure(c(list(target = "difference"), specification), class = "flexbayes_prior_update"))
  }
  if (has_outcome) {
    if (!missing(remove)) stop("remove applies only to parameter beliefs; outcome removal is not supported.")
    specification = prior_outcome_update(outcome, beliefs)
    return(structure(c(list(target = "outcome"), specification), class = "flexbayes_prior_update"))
  }

  if (!is.character(parameter) || length(parameter) != 1L || is.na(parameter) || !nzchar(parameter)) stop("parameter must be one nonempty parameter name.")
  prior_validate_removals(remove)
  if (!length(beliefs) && !length(remove)) stop("Supply at least one belief() object or removal index.")
  if (any(!vapply(beliefs, inherits, logical(1), "prior_belief"))) stop("All additions must be belief objects.")
  structure(list(target = "parameter", parameter = parameter, beliefs = beliefs, remove = remove), class = "flexbayes_prior_update")
}

#' @export
`+.flexbayes_priors` = function(e1, e2) {
  if (missing(e2)) stop("Use priors + update_priors(...); unary + is not supported.")
  if (!inherits(e1, "flexbayes_priors") || !inherits(e2, "flexbayes_prior_update")) {
    stop("Use a flexbayes_priors object on the left and an update_priors() specification on the right.")
  }

  if (e2$target %in% c("outcome", "difference")) return(prior_maybe_switch_package(prior_apply_outcome(e1, e2)))
  if (!identical(e2$target, "parameter")) stop("Unknown update target; construct a specification with update_priors().")

  # Unname the list so a belief argument cannot masquerade as parameter/remove.
  remove = e2$remove
  object = do.call(prior_apply_parameter, c(list(object = e1, parameter = e2$parameter, remove = remove), unname(e2$beliefs)))
  prior_maybe_switch_package(object)
}

# Use the identical method for right-hand dispatch as well, so unsupported
# operands receive the same error and our two classes cannot conflict in Ops.
#' @export
`+.flexbayes_prior_update` = `+.flexbayes_priors`

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
