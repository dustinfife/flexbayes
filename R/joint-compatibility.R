# Conservative compatibility checks. Bounds describe an outer set of feasible
# Normal moments, so an empty outer set proves a conflict; a nonempty set does
# not establish feasibility of the complete mean/variance problem.
prior_joint_moment_bounds = function(problem) {
  scale = sqrt(as.vector(problem$X^2 %*% problem$sd^2))
  reference = as.vector(problem$X %*% problem$mean)
  lower = rep(0,nrow(problem$X)); upper = rep(Inf,nrow(problem$X))
  sources = vector("list",nrow(problem$X))
  for (i in seq_len(nrow(problem$X))) {
    same = which(apply(problem$X,1,function(row) identical(unname(row),unname(problem$X[i,]))))
    indices = which(vapply(problem$outcome,function(r) r$configuration %in% same,logical(1)))
    lo = hi = matrix(numeric(),ncol=3,dimnames=list(NULL,c("Intercept","SD_slope","Record")))
    for (index in indices) {
      r = problem$outcome[[index]]
      if (r$type %in% c("mean","mode")) {
        lo = rbind(lo,c((r$target-r$tolerance-reference[i])/scale[i],0,index))
        hi = rbind(hi,c((r$target+r$tolerance-reference[i])/scale[i],0,index))
        next
      }
      h = (r$threshold-reference[i])/scale[i]
      q = stats::qnorm(c(max(0,r$target-r$tolerance),min(1,r$target+r$tolerance)))
      if (r$operator %in% c("<","<=")) {
        if (is.finite(q[2])) lo = rbind(lo,c(h,-q[2],index))
        if (is.finite(q[1])) hi = rbind(hi,c(h,-q[1],index))
      } else {
        if (is.finite(q[1])) lo = rbind(lo,c(h,q[1],index))
        if (is.finite(q[2])) hi = rbind(hi,c(h,q[2],index))
      }
    }
    for (a in seq_len(nrow(lo))) for (b in seq_len(nrow(hi))) {
      slope = lo[a,2]-hi[b,2]
      rhs = hi[b,1]-lo[a,1]
      if (slope==0 && rhs < -1e-12) return(list(conflict=same,records=as.integer(c(lo[a,3],hi[b,3])),lines=list(lower=lo,upper=hi)))
      if (slope==0) next
      bound = rhs/slope
      # Widen, never tighten, the derived outer set for numerical roundoff.
      margin = 1e-12*max(1,abs(bound))
      if (slope>0 && bound+margin<upper[i]) {
        upper[i] = bound+margin
        sources[[i]] = union(sources[[i]],as.integer(c(lo[a,3],hi[b,3])))
      }
      if (slope<0 && bound-margin>lower[i]) {
        lower[i] = bound-margin
        sources[[i]] = union(sources[[i]],as.integer(c(lo[a,3],hi[b,3])))
      }
    }
    if (upper[i]<lower[i]) return(list(conflict=same,records=indices,lower_sd=lower[i]*scale[i],upper_sd=upper[i]*scale[i],lines=list(lower=lo,upper=hi)))
  }
  list(lower=lower,upper=upper,scale=scale,records=sources)
}

prior_joint_variance_distance = function(v, A, variance_lower, variance_upper) {
  values = as.vector(A%*%v)
  residual = values-pmin(variance_upper,pmax(variance_lower,values))
  list(value=sum(residual^2)/2,gradient=as.vector(crossprod(A,residual)),residual=residual)
}
prior_joint_variance_value = function(v, ...) prior_joint_variance_distance(v, ...)$value
prior_joint_variance_gradient = function(v, ...) prior_joint_variance_distance(v, ...)$gradient

# A verified separating hyperplane: A' r >= 0, but sup_{y in box} r'y < 0.
# These inequalities contradict A v in box for every v >= 0. Optimization only
# searches for a witness; its convergence flag is never an impossibility proof.
prior_joint_variance_certificate = function(A, lower, upper) {
  fit = tryCatch(stats::optim(rep(1,ncol(A)),prior_joint_variance_value,
    prior_joint_variance_gradient,A=A,variance_lower=lower,variance_upper=upper,method="L-BFGS-B",
    lower=rep(0,ncol(A)),control=list(maxit=2000,factr=10,pgtol=1e-12)),error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  r = prior_joint_variance_distance(fit$par,A,lower,upper)$residual
  finite = which(is.finite(upper))
  if (!length(finite) || !any(r!=0)) return(NULL)
  support = colSums(A[finite,,drop=FALSE])
  dot = as.vector(crossprod(A,r))
  epsilon = 128*.Machine$double.eps*(nrow(A)+ncol(A)+1)
  # A tiny correction supported only by finite upper bounds makes the dual
  # inequalities conservative without introducing infinite box products.
  correction = if (any(support>0)) max(0,(epsilon-dot[support>0])/support[support>0]) else 0
  r[finite] = r[finite]+correction
  dot = as.vector(crossprod(A,r))
  error = epsilon*colSums(abs(A*r))
  if (any(dot<error)) return(NULL)
  positive = which(r>0); negative = which(r<0)
  if (any(!is.finite(upper[positive]))) return(NULL)
  terms = c(r[positive]*upper[positive],r[negative]*lower[negative])
  gap = -sum(terms)
  if (!is.finite(gap) || gap<=epsilon*max(1,sum(abs(terms)))) return(NULL)
  list(weights=r,nonnegative=dot,gap=gap,matrix=A,lower=lower,upper=upper)
}

prior_joint_compatibility = function(problem) {
  basic = prior_joint_basic_compatibility(problem)
  if (!is.null(basic)) return(basic)
  bounds = prior_joint_moment_bounds(problem)
  if (!is.null(bounds$conflict)) return(list(issue="location and variance",
    message="Outcome beliefs are incompatible with a Normal expected-outcome distribution at this configuration.",
    configurations=bounds$conflict,records=bounds$records,
    explanation="These answers cannot all describe the location and uncertainty of one Normal expected outcome, even within their acceptance intervals.",
    bounds=if (is.null(bounds$lower_sd)) NULL else data.frame(Configuration=bounds$conflict,
      Lower_SD=bounds$lower_sd,Upper_SD=bounds$upper_sd),
    technical=c(list(method="Normal moment half-space bounds"),bounds$lines)))
  squared = problem$X^2
  for (i in seq_len(nrow(squared))) for (j in seq_len(nrow(squared))) {
    if (i==j || any(squared[j,]==0 & squared[i,]>0)) next
    active = squared[j,]>0
    cap = max(squared[i,active]/squared[j,active])
    required = bounds$lower[i]*bounds$scale[i]
    allowed = sqrt(cap)*bounds$upper[j]*bounds$scale[j]
    if (!is.finite(allowed) || required<=allowed*(1+1e-10)) next
    records = unique(c(bounds$records[[i]],bounds$records[[j]]))
    return(list(issue="variance",message="Outcome beliefs are incompatible: the requested uncertainty ratio exceeds this model's independent-Normal variance restriction.",
      configurations=c(j,i),records=records,
      ratio=list(numerator=i,denominator=j,required=required/(bounds$upper[j]*bounds$scale[j]),maximum=sqrt(cap)),
      bounds=data.frame(Configuration=c(j,i),Lower_SD=bounds$lower[c(j,i)]*bounds$scale[c(j,i)],Upper_SD=bounds$upper[c(j,i)]*bounds$scale[c(j,i)]),
      technical=list(method="Nonnegative variance ratio bound",squared_design=squared,variance_ratio_limit=cap)))
  }
  A = sweep(squared,2,problem$sd^2,"*")/bounds$scale^2
  proof = prior_joint_variance_certificate(A,bounds$lower^2,bounds$upper^2)
  if (is.null(proof)) return(NULL)
  configurations = which(abs(proof$weights)>0)
  list(issue="variance",message="Outcome beliefs are incompatible: no nonnegative coefficient variances can produce all the requested expected-outcome uncertainties.",
    configurations=configurations,records=unique(unlist(bounds$records[configurations])),
    explanation="The requested uncertainties cannot be combined across these configurations using independent Normal coefficient priors, even within the acceptance intervals.",
    bounds=data.frame(Configuration=configurations,Lower_SD=bounds$lower[configurations]*bounds$scale[configurations],Upper_SD=bounds$upper[configurations]*bounds$scale[configurations]),
    technical=c(list(method="Separating hyperplane for nonnegative coefficient variances"),proof))
}

# Preserve the string interface used by existing internal callers and tests.
prior_joint_conflict = function(problem) {
  certificate = prior_joint_compatibility(problem)
  if (is.null(certificate)) return(NULL)
  certificate$message
}

prior_joint_summary_spread_bounds = function(bounds) {
  for (i in seq_len(nrow(bounds))) {
    lower = prior_joint_summary_number(bounds$Lower_SD[i])
    upper = prior_joint_summary_number(bounds$Upper_SD[i])
    text = if (!is.finite(bounds$Upper_SD[i])) paste("at least",lower) else if (lower==upper)
      paste("about",lower) else paste("between",lower,"and",upper)
    if (bounds$Lower_SD[i]>bounds$Upper_SD[i]) text = paste("at least",lower,"but at most",upper)
    prior_console_line(paste0("Implied uncertainty (SD) at configuration ",bounds$Configuration[i],": ",text,"."))
  }
}

prior_joint_summary_conflict = function(object, detailed=FALSE) {
  certificate = object$joint$compatibility
  if (is.null(certificate)) return(invisible(NULL))
  prior_console_line(paste("Conflict concerns:",certificate$issue))
  for (configuration in certificate$configurations) {
    records = which(vapply(object$joint$problem$outcome,function(r) r$configuration==configuration,logical(1)))
    numbers = match(intersect(certificate$records,records),records)
    specification = object$outcomes[[configuration]]
    at = if (length(specification$at)) paste(names(specification$at),
      vapply(specification$at,prior_joint_summary_number,character(1)),sep=" = ",collapse=", ") else "intercept only"
    prior_console_line(paste0("Configuration ",configuration," (",at,"): beliefs ",paste(numbers,collapse=", "),"."))
  }
  if (!is.null(certificate$ratio)) {
    ratio = certificate$ratio
    required = format(signif(ratio$required,6),trim=TRUE)
    maximum = format(signif(ratio$maximum,6),trim=TRUE)
    prior_console_line(paste0("These beliefs require at least ",required,
      " times as much uncertainty (SD) at configuration ",ratio$numerator,
      " as at configuration ",ratio$denominator,". This model allows at most ",maximum," times."))
    prior_console_line("The restriction comes from independent Normal coefficient priors and the model's original centering values.")
  } else {
    if (!is.null(certificate$bounds)) prior_joint_summary_spread_bounds(certificate$bounds)
    if (!is.null(certificate$explanation)) prior_console_line(certificate$explanation)
  }
  prior_console_line("Your answers are retained; no verified adjusted priors were selected.")
  if (detailed) {
    if (!is.null(certificate$bounds)) print(certificate$bounds,row.names=FALSE,digits=8)
    cat("Mathematical certificate:\n")
    print(certificate$technical)
  }
  invisible(NULL)
}
