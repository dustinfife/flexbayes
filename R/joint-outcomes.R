#' Elicit the expected outcome at a predictor configuration
#'
#' Add this specification to an identify_priors() object with `+`. Beliefs concern
#' the conditional expected response, excluding residual noise. The original
#' object and coefficient belief statements are never modified in place.
#' @param at Named list specifying one finite numeric value for every predictor.
#'   For an intercept-only model use list().
#' @param ... Mean, mode or probability inequality objects from belief(), using
#'   scale="original". Multiple specifications are reconciled jointly.
#' @details Initially supports Gaussian identity models with an intercept,
#' numeric predictors, and selected independent Normal coefficient priors.
#' Student's t priors, categorical predictors, offsets and residual-scale
#' elicitation are unsupported. The existing centered-intercept parameterization
#' is retained; independence refers to that parameterization.
#'
#' Reconciliation has three priorities: satisfy outcome beliefs within their
#' tolerances; minimize the largest coefficient-belief violation beyond its
#' tolerance; then minimize KL(adjusted prior || baseline prior). Defaults enter
#' only the last tier. This is a deterministic constrained numerical search,
#' not a guarantee of global optimality. An unverified solution is never exported.
#' Original coefficient-only selections define the baseline on every update,
#' preventing successive outcome additions from moving the baseline.
#'
#' The fixed-variance Mahalanobis projection initializes the numerical search.
#' General reconciliation may change both means and SDs. Underidentified problems
#' are resolved by the stated preference criterion, not by claiming identification.
#' Relaxed coefficient beliefs produce a warning and remain visible in
#' summary(..., diagnostics=TRUE). Incompatible or numerically unresolved outcome
#' specifications are retained but block check_prior(), as_rstanarm() and fit().
#' @return A flexbayes_outcome_update specification; addition returns revised priors.
#' @examples
#' change = update_outcome(at=list(x=2), belief(mode=10),
#'   belief(value < 8, probability=.05))
#' \dontrun{
#' d = data.frame(x=1:6, y=c(3, 5, 4, 8, 9, 11))
#' p = identify_priors(y ~ x, d) + change
#' summary(p, diagnostics=TRUE)
#' }
#' @export
update_outcome = function(at, ...) {
  if (!is.list(at) || length(at) && (is.null(names(at)) || anyNA(names(at)) ||
      any(!nzchar(names(at))) || anyDuplicated(names(at)))) stop("at must be a uniquely named list.")
  if (any(!vapply(at, function(x) is.numeric(x) && !is.complex(x) && length(x)==1L &&
      is.null(dim(x)) && is.finite(x), logical(1)))) stop("Each predictor configuration value must be one finite number.")
  records = list(...)
  if (!length(records) || any(!vapply(records, inherits, logical(1), "prior_belief"))) stop("Supply belief() objects.")
  if (any(!vapply(records, function(r) r$scale=="original" && r$type %in% c("mean", "mode", "probability"), logical(1))))
    stop("Expected-outcome elicitation supports original-scale mean, mode and probability beliefs only.")
  structure(list(at=at, beliefs=unname(records)), class="flexbayes_outcome_update")
}

prior_apply_outcome = function(object, specification) {
  for (i in seq_along(specification$beliefs)) {
    specification$beliefs[[i]]$id = object$next_id
    object$next_id = object$next_id + 1L
  }
  object$outcomes = c(object$outcomes, list(unclass(specification)))
  prior_select_compatible(object)
}

prior_joint_design = function(object, specifications) {
  if (object$model$family$family != "gaussian" || object$model$family$link != "identity")
    stop("Joint outcome elicitation currently supports Gaussian identity models only.")
  if (!"(Intercept)" %in% colnames(object$model$X)) stop("Joint outcome elicitation requires an intercept.")
  variables = all.vars(stats::delete.response(object$model$terms))
  response = attr(object$model$terms,"response")
  predictors = if (response) object$model$frame[-response] else object$model$frame
  if (any(!vapply(object$model$data[variables],is.numeric,logical(1))) ||
      any(!vapply(predictors,is.numeric,logical(1))))
    stop("Joint outcome elicitation currently requires numeric predictors, including their formula transformations.")
  rows = lapply(specifications, function(specification) {
    if (!setequal(names(specification$at), variables)) stop("at must specify exactly these predictors: ", paste(variables, collapse=", "))
    data = if (length(variables)) as.data.frame(specification$at) else data.frame(row.names=1L)
    X = prior_prediction_matrix(object, data)
    slopes = setdiff(colnames(X), "(Intercept)")
    X[,slopes] = sweep(X[,slopes,drop=FALSE], 2, object$model$reference[slopes], "-")
    X
  })
  do.call(rbind, rows)
}

# Parameters are dimensionless mean displacements and log-SD ratios.
prior_joint_parameters = function(theta, problem) {
  n = length(problem$mean)
  list(mean=problem$mean + problem$sd*theta[seq_len(n)],
    sd=problem$sd*exp(theta[n+seq_len(n)]))
}

prior_joint_values = function(theta, problem, coefficient=FALSE) {
  parameters = prior_joint_parameters(theta, problem)
  records = if (coefficient) problem$coefficient else problem$outcome
  p = length(parameters$mean)
  value = numeric(length(records))
  gradient = matrix(0,length(records),2*p)
  for (i in seq_along(records)) {
    record = records[[i]]
    row = if (coefficient) diag(p)[record$parameter,] else problem$X[record$configuration,]
    mean = sum(row*parameters$mean)
    variance = sum((row*parameters$sd)^2)
    sd = sqrt(variance)
    mean_gradient = row*problem$sd
    if (record$type %in% c("mean","mode")) {
      value[i] = mean
      gradient[i,seq_len(p)] = mean_gradient
    } else if (record$type == "sd") {
      value[i] = sd
      gradient[i,p+seq_len(p)] = (row*parameters$sd)^2/sd
    } else {
      z = (record$threshold-mean)/sd
      sign = if (record$operator %in% c("<","<=")) 1 else -1
      value[i] = if (sign==1) stats::pnorm(z) else stats::pnorm(z,lower.tail=FALSE)
      gradient[i,seq_len(p)] = -sign*stats::dnorm(z)*mean_gradient/sd
      gradient[i,p+seq_len(p)] = -sign*stats::dnorm(z)*z*(row*parameters$sd)^2/variance
    }
  }
  list(value=value, gradient=gradient)
}

prior_joint_kl = function(theta, problem) {
  p = length(problem$mean)
  u = theta[seq_len(p)]; l = theta[p+seq_len(p)]
  list(value=sum(u^2+expm1(2*l)-2*l)/2, gradient=c(u,expm1(2*l)))
}

# Scalar objective and inequalities (g <= 0), with analytic derivatives.
prior_joint_system = function(theta, problem, tier, coefficient_limit=0, coefficient_scale=1) {
  p = length(problem$mean)*2L
  parameters = theta[seq_len(p)]
  outcome = prior_joint_values(parameters, problem)
  target = vapply(problem$outcome, `[[`, numeric(1), "target")
  tolerance = vapply(problem$outcome, `[[`, numeric(1), "tolerance")
  difference = outcome$value-target
  # A small interior margin protects final checks from floating-point rounding.
  # Retry without it when the tolerance boxes have only a narrow intersection.
  width = tolerance*(1-problem$interior)
  g = c((difference-width)/problem$units,(-difference-width)/problem$units)
  jacobian = rbind(outcome$gradient/problem$units,-outcome$gradient/problem$units)
  coefficient = prior_joint_values(parameters, problem, TRUE)
  if (length(problem$coefficient) && tier != 1L) {
    ct = vapply(problem$coefficient, `[[`, numeric(1), "target")
    tol = vapply(problem$coefficient, `[[`, numeric(1), "tolerance")
    delta = coefficient$value-ct
    allowance = tol*(1+if (tier==2L) coefficient_scale*theta[p+1L] else coefficient_limit)
    g = c(g,(delta-allowance)/problem$coefficient_units,(-delta-allowance)/problem$coefficient_units)
    jacobian = rbind(jacobian,coefficient$gradient/problem$coefficient_units,-coefficient$gradient/problem$coefficient_units)
  }
  if (tier==2L) {
    jacobian = cbind(jacobian, c(rep(0,2*length(target)),rep(-tol*coefficient_scale/problem$coefficient_units,2)))
    g = c(g,-theta[p+1L])
    jacobian = rbind(jacobian,c(rep(0,p),-1))
    return(list(value=theta[p+1L],gradient=c(rep(0,p),1),g=g,jacobian=jacobian))
  }
  kl = prior_joint_kl(parameters,problem)
  list(value=kl$value,gradient=kl$gradient,g=g,jacobian=jacobian)
}

prior_joint_augmented = function(theta, problem, tier, coefficient_limit, coefficient_scale, lambda, rho) {
  system = prior_joint_system(theta,problem,tier,coefficient_limit,coefficient_scale)
  if (any(!is.finite(c(system$value,system$gradient,system$g,system$jacobian))))
    return(list(value=1e100,gradient=rep(0,length(theta))))
  multiplier = pmax(0,lambda+rho*system$g)
  list(value=system$value+sum(multiplier^2-lambda^2)/(2*rho),
    gradient=system$gradient+as.vector(crossprod(system$jacobian,multiplier)))
}
prior_joint_augmented_value = function(theta, ...) prior_joint_augmented(theta, ...)$value
prior_joint_augmented_gradient = function(theta, ...) prior_joint_augmented(theta, ...)$gradient

prior_joint_solve = function(start, problem, tier=3L, coefficient_limit=0, coefficient_scale=1) {
  lambda = rep(0,length(prior_joint_system(start,problem,tier,coefficient_limit,coefficient_scale)$g))
  rho = 10
  previous = Inf
  for (iteration in seq_len(24L)) {
    fit = stats::optim(start,prior_joint_augmented_value,prior_joint_augmented_gradient,
      problem=problem,tier=tier,coefficient_limit=coefficient_limit,coefficient_scale=coefficient_scale,
      lambda=lambda,rho=rho,method="BFGS",control=list(maxit=600,reltol=1e-13))
    current = prior_joint_system(fit$par,problem,tier,coefficient_limit,coefficient_scale)
    violation = max(c(0,current$g))
    step = max(abs(fit$par-start))
    lambda = pmax(0,lambda+rho*current$g)
    start = fit$par
    if (violation < 1e-12 && step < 1e-7 && fit$convergence==0L) break
    if (violation > previous*.3) rho = min(rho*10,1e14)
    previous = violation
  }
  # Restore numerical feasibility after the augmented-Lagrangian iterations.
  # Baseline-scaled constraints can be tiny in dollar-valued examples; an
  # optimizer's absolute stopping rule is insufficient for belief tolerances.
  for (polish in seq_len(20L)) {
    current = prior_joint_system(start,problem,tier,coefficient_limit,coefficient_scale)
    active = which(current$g > 0)
    if (!length(active)) break
    J = current$jacobian[active,,drop=FALSE]
    singular = svd(J)
    keep = singular$d > max(singular$d)*1e-12
    if (!any(keep)) break
    correction = as.vector(singular$v[,keep,drop=FALSE] %*%
      (crossprod(singular$u[,keep,drop=FALSE],current$g[active])/singular$d[keep]))
    start = start-correction
  }
  current = prior_joint_system(start,problem,tier,coefficient_limit,coefficient_scale)
  list(par=start,convergence=fit$convergence,iterations=iteration,
    violation=max(c(0,current$g)),objective=current$value,
    stationarity=prior_joint_stationarity(current))
}

# A nonnegative multiplier certificate checks local first-order optimality.
# Dropping negative multipliers is conservative: failure leaves the problem
# unresolved rather than accepting a merely feasible, nonstationary solution.
prior_joint_stationarity = function(system) {
  active = which(system$g >= -1e-7)
  residual = system$gradient
  while (length(active)) {
    J = system$jacobian[active,,drop=FALSE]
    singular = svd(t(J))
    keep = singular$d > max(singular$d)*1e-12
    if (!any(keep)) break
    multiplier = as.vector(singular$v[,keep,drop=FALSE] %*%
      (crossprod(singular$u[,keep,drop=FALSE],-system$gradient)/singular$d[keep]))
    if (all(multiplier>=0)) {
      residual = c(system$gradient+as.vector(crossprod(J,multiplier)),
        multiplier*system$g[active])
      break
    }
    active = active[-which.min(multiplier)]
  }
  max(abs(residual))
}

prior_joint_projection = function(X, target, mean, sd) {
  A = sweep(X,2,sd,"*")
  singular = svd(A)
  keep = singular$d > max(singular$d)*1e-12
  adjustment = if (any(keep)) as.vector(singular$v[,keep,drop=FALSE] %*%
    (crossprod(singular$u[,keep,drop=FALSE],target-as.vector(X%*%mean))/singular$d[keep])) else rep(0,length(mean))
  mean+sd*adjustment
}

prior_joint_problem = function(object) {
  X = prior_joint_design(object,object$outcomes)
  names = colnames(X)
  baseline = object$selected
  if (any(vapply(baseline[names],is.null,logical(1)))) stop("Resolve coefficient-only priors before adding outcome beliefs.")
  if (any(vapply(baseline[names],function(p) p$distribution != "normal",logical(1))))
    stop("Joint outcome elicitation requires independent Normal priors, including the centered intercept; Student's t priors are unsupported.")
  mean = vapply(baseline[names],`[[`,numeric(1),"mean")
  sd = vapply(baseline[names],`[[`,numeric(1),"sd")
  outcome = list(); coefficient = list()
  for (i in seq_along(object$outcomes)) for (record in object$outcomes[[i]]$beliefs) {
    record$configuration = i
    outcome = c(outcome,list(record))
  }
  for (i in seq_along(names)) for (record in object$beliefs[[names[i]]]) {
    record$parameter = i
    coefficient = c(coefficient,list(record))
  }
  variance = as.vector(X^2 %*% sd^2)
  units = vapply(outcome,function(r) if (r$type=="probability") 1 else sqrt(variance[r$configuration]),numeric(1))
  coefficient_units = vapply(coefficient,function(r) if (r$type=="probability") 1 else sd[r$parameter],numeric(1))
  list(X=X,mean=mean,sd=sd,outcome=outcome,coefficient=coefficient,
    units=units,coefficient_units=coefficient_units,baseline=baseline,interior=.001)
}

prior_joint_conflict = function(problem) {
  for (i in seq_len(nrow(problem$X))) {
    rows = which(apply(problem$X,1,function(x) identical(unname(x),unname(problem$X[i,]))))
    records = Filter(function(r) r$configuration %in% rows,problem$outcome)
    check = prior_check_probability_constraints(records)
    if (check$contradictory) return(check$message)
  }
  location = Filter(function(r) r$type %in% c("mean","mode"),problem$outcome)
  if (!length(location)) return(NULL)
  A = problem$X[vapply(location,`[[`,integer(1),"configuration"),,drop=FALSE]
  b = vapply(location,`[[`,numeric(1),"target")
  tol = vapply(location,`[[`,numeric(1),"tolerance")
  projected = prior_joint_projection(A,b,problem$mean,problem$sd)
  residual = b-as.vector(A%*%projected)
  # A left-null-space certificate proves that the tolerance boxes cannot meet.
  if (max(abs(crossprod(A,residual))) < 1e-8*max(1,sum(abs(A)*max(abs(residual)))) &&
      abs(sum(residual*b)) > sum(abs(residual)*tol)+1e-9*max(1,sum(abs(residual*b))))
    return("Outcome location constraints are incompatible with the model matrix, even within their tolerances.")
  NULL
}

prior_joint_details = function(problem,theta,coefficient=FALSE) {
  records = if (coefficient) problem$coefficient else problem$outcome
  values = prior_joint_values(theta,problem,coefficient)$value
  target = vapply(records,`[[`,numeric(1),"target")
  tolerance = vapply(records,`[[`,numeric(1),"tolerance")
  residual = (values-target)/tolerance
  data.frame(id=vapply(records,`[[`,numeric(1),"id"),target=target,achieved=values,
    tolerance=tolerance,residual=residual,satisfied=is.finite(residual) & abs(residual)<=1)
}

prior_joint_reconcile = function(object) {
  object$joint = NULL
  problem = prior_joint_problem(object)
  conflict = prior_joint_conflict(problem)
  base = rep(0,2*length(problem$mean))
  joint = list(baseline=problem$baseline,beliefs=object$beliefs,
    specifications=object$outcomes,design=problem$X,problem=problem,
    status="unresolved",message=conflict)
  if (!is.null(conflict)) {
    joint$status = "incompatible"
    object$joint = joint
    warning(conflict,call.=FALSE)
    return(object)
  }
  # A projection provides a deterministic, unit-aware starting point, including
  # single probability constraints that do not identify coefficient parameters.
  targets = vapply(problem$outcome,function(r) {
    if (r$type != "probability") return(r$target)
    p = if (r$operator %in% c("<","<=")) r$target else 1-r$target
    if (p<=0 || p>=1) stop("Joint Normal outcome beliefs require probabilities strictly between zero and one.")
    r$threshold-stats::qnorm(p)*sqrt(sum((problem$X[r$configuration,]*problem$sd)^2))
  },numeric(1))
  A = problem$X[vapply(problem$outcome,`[[`,integer(1),"configuration"),,drop=FALSE]
  location = prior_joint_projection(A,targets,problem$mean,problem$sd)
  start = c((location-problem$mean)/problem$sd,rep(0,length(problem$sd)))
  initial = prior_joint_solve(start,problem,tier=1L)
  if (!all(prior_joint_details(problem,initial$par)$satisfied) || initial$violation>1e-8) {
    problem$interior = 0
    initial = prior_joint_solve(start,problem,tier=1L)
  }
  joint$problem = problem
  stages = list(outcome=initial)
  theta = initial$par
  outcome = prior_joint_details(problem,theta)
  if (all(outcome$satisfied) && length(problem$coefficient)) {
    coefficient = prior_joint_details(problem,theta,TRUE)
    scale = max(1,pmax(abs(coefficient$residual)-1,0))
    second = prior_joint_solve(c(theta,1),problem,tier=2L,coefficient_scale=scale)
    theta2 = head(second$par,length(base))
    details = prior_joint_details(problem,theta2,TRUE)
    limit = max(c(0,abs(details$residual)-1))
    third = prior_joint_solve(theta2,problem,tier=3L,coefficient_limit=limit*(1+1e-8))
    stages = c(stages,list(coefficient=second,KL=third))
    theta = third$par
  }
  joint$outcome = prior_joint_details(problem,theta)
  joint$coefficient = prior_joint_details(problem,theta,TRUE)
  joint$stages = stages
  joint$theta = theta
  joint$kl = prior_joint_kl(theta,problem)$value
  joint$reference_used = TRUE
  joint$underdetermined = qr(prior_joint_values(theta,problem)$gradient/problem$units)$rank < length(theta)
  verified = all(joint$outcome$satisfied) && all(vapply(stages,function(s) s$convergence==0 && is.finite(s$objective) &&
    s$violation<1e-8 && is.finite(s$stationarity) && s$stationarity<1e-4,logical(1)))
  if (!verified) {
    joint$message = "No converged feasible joint prior was verified. Outcome beliefs may conflict, or numerical reconciliation failed; revise the beliefs or tolerances."
    object$joint = joint
    warning(joint$message,call.=FALSE)
    return(object)
  }
  parameters = prior_joint_parameters(theta,problem)
  for (i in seq_along(problem$mean)) object$selected[[names(problem$mean)[i]]] =
    list(distribution="normal",mean=unname(parameters$mean[i]),sd=unname(parameters$sd[i]))
  joint$adjusted = object$selected
  joint$status = if (all(joint$coefficient$satisfied)) "compatible" else "reconciled"
  object$joint = joint
  if (joint$status=="reconciled") warning("Some coefficient beliefs were relaxed to satisfy outcome beliefs; inspect summary(..., diagnostics=TRUE).",call.=FALSE)
  object
}

# Recompute achieved quantities; do not trust stale reconciliation annotations.
prior_joint_validate = function(object) {
  joint = object$joint
  if (is.null(joint)) return(invisible(object))
  if (!joint$status %in% c("compatible","reconciled")) stop("Outcome reconciliation is ",joint$status,": ",joint$message)
  if (!identical(object$beliefs,joint$beliefs) || !identical(object$outcomes,joint$specifications) ||
      !isTRUE(all.equal(object$selected,joint$adjusted))) stop("Joint prior metadata or selected distributions changed; re-elicit the priors.")
  prior_validate_standardization_data(object)
  baseline = object
  baseline$selected = joint$baseline
  rebuilt = prior_joint_problem(baseline)
  fields = setdiff(names(rebuilt),"interior")
  if (!isTRUE(all.equal(rebuilt[fields],joint$problem[fields])))
    stop("Reconciliation baseline or constraints changed; re-elicit the priors.")
  parameters = prior_joint_parameters(joint$theta,joint$problem)
  for (i in seq_along(parameters$mean)) {
    selected = object$selected[[names(parameters$mean)[i]]]
    if (!isTRUE(all.equal(unname(c(selected$mean,selected$sd)),
        unname(c(parameters$mean[i],parameters$sd[i])))))
      stop("Adjusted priors do not match the reconciliation parameters.")
  }
  X = prior_prediction_matrix(object,object$model$data)
  if (!isTRUE(all.equal(X,object$model$X))) stop("Model data changed since outcome elicitation; identify priors again.")
  if (!isTRUE(all.equal(prior_joint_design(object,object$outcomes),joint$design))) stop("Outcome prediction design changed; re-elicit the priors.")
  outcome = prior_joint_details(joint$problem,joint$theta)
  coefficient = prior_joint_details(joint$problem,joint$theta,TRUE)
  if (!all(outcome$satisfied) || !isTRUE(all.equal(outcome,joint$outcome)) ||
      !isTRUE(all.equal(coefficient,joint$coefficient))) stop("Reconciliation constraints or diagnostics do not validate.")
  invisible(object)
}

prior_joint_summary = function(x) {
  object = x$object; joint = object$joint
  print(object)
  cat("\nOutcome reconciliation:",joint$status,"\n")
  cat("Beliefs concern conditional expected outcomes, without residual noise.\n")
  if (!is.null(joint$message)) prior_console_line(joint$message)
  if (!x$diagnostics) return(invisible(x))
  index = 0L
  for (specification in object$outcomes) {
    cat("\nAt:",if (length(specification$at)) paste(names(specification$at),unlist(specification$at),sep="=",collapse=", ") else "intercept only","\n")
    for (record in specification$beliefs) {
      index = index+1L
      prior_console_line(record$label)
      if (!is.null(joint$outcome)) print(joint$outcome[index,,drop=FALSE],row.names=FALSE,digits=8)
    }
  }
  cat("\nCoefficient reconciliation (centered-intercept parameterization):\n")
  for (name in names(joint$problem$mean)) {
    before = joint$baseline[[name]]; after = object$selected[[name]]
    cat(name,":",prior_format(before),"->",prior_format(after),"\n")
    cat("  Change in mean:",after$mean-before$mean,"; change in SD:",after$sd-before$sd,"\n")
    records = object$beliefs[[name]]
    for (record in records) {
      row = if (is.null(joint$coefficient)) data.frame() else
        joint$coefficient[joint$coefficient$id==record$id,,drop=FALSE]
      prior_console_line(paste(record$label,if (nrow(row) && row$satisfied) "[satisfied]" else "[relaxed / unresolved]"))
      if (nrow(row)) print(row,row.names=FALSE,digits=8)
    }
  }
  if (!is.null(joint$kl)) cat("\nKL(adjusted || baseline):",joint$kl,"\n")
  cat("Reference anchoring: coefficient-only selections, with backend defaults where not elicited.\n")
  cat("Underdetermined constraints use the stated minimax/KL criterion; numerical search does not certify a unique global optimum.\n")
  if (length(joint$stages)) print(do.call(rbind,lapply(joint$stages,function(s) data.frame(
    convergence=s$convergence,iterations=s$iterations,scaled_violation=s$violation,
    stationarity=s$stationarity,objective=s$objective))))
  invisible(x)
}
