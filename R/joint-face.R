# The fixed minimax face has linear coefficient constraints in means and SDs.
# Probability bands are equivalent Normal-quantile inequalities, not new beliefs.
prior_joint_face_add = function(system, values, jacobian, curvature = NULL) {
  system$g = c(system$g, values)
  system$jacobian = rbind(system$jacobian, jacobian)
  if (is.null(curvature)) curvature = rep(list(matrix(0, ncol(jacobian), ncol(jacobian))), length(values))
  system$curvature = c(system$curvature, curvature)
  system
}

prior_joint_face_system = function(x, face) {
  problem = face$problem
  n = length(problem$mean)
  indices = seq_len(n)
  r = x[n + indices]
  mean = problem$mean + problem$sd * x[indices]
  sd = face$sd * r
  log_ratio = log(face$sd / problem$sd) + log(r)
  system = list(value = sum(x[indices]^2 + expm1(2 * log_ratio) - 2 * log_ratio) / 2,
    gradient = c(x[indices], (face$sd / problem$sd)^2 * r - 1 / r),
    hessian = diag(c(rep(1, n), (face$sd / problem$sd)^2 + 1 / r^2)),
    g = numeric(), jacobian = matrix(numeric(), 0, 2 * n), curvature = list())

  for (record in problem$outcome) {
    row = problem$X[record$configuration, ]
    location = sum(row * mean)
    variance_weights = (row * face$sd)^2
    spread = sqrt(sum(variance_weights * r^2))
    unit = sqrt(sum((row * problem$sd)^2))
    dm = row * problem$sd
    ds = variance_weights * r / spread
    width = record$tolerance * (1 - problem$interior)
    if (record$type %in% c("mean", "mode")) {
      system = prior_joint_face_add(system,
        c(location - record$target - width, record$target - location - width) / unit,
        rbind(c(dm, rep(0, n)), c(-dm, rep(0, n))) / unit)
      next
    }
    sign = if (record$operator %in% c("<", "<=")) 1 else -1
    quantiles = stats::qnorm(c(max(0, record$target - width), min(1, record$target + width)))
    H = matrix(0, 2 * n, 2 * n)
    H[n + indices, n + indices] = (diag(variance_weights, nrow = n) / spread -
      tcrossprod(variance_weights * r) / spread^3) / unit
    if (is.finite(quantiles[2])) system = prior_joint_face_add(system,
      (sign * (record$threshold - location) - quantiles[2] * spread) / unit,
      matrix(c(-sign * dm, -quantiles[2] * ds) / unit, 1), list(-quantiles[2] * H))
    if (is.finite(quantiles[1])) system = prior_joint_face_add(system,
      (quantiles[1] * spread - sign * (record$threshold - location)) / unit,
      matrix(c(sign * dm, quantiles[1] * ds) / unit, 1), list(quantiles[1] * H))
  }

  system$coefficient_start = length(system$g)
  for (record in problem$coefficient) {
    j = record$parameter
    dm = ds = numeric(n)
    dm[j] = problem$sd[j]
    ds[j] = face$sd[j]
    unit = problem$sd[j]
    if (record$type != "probability") {
      value = if (record$type == "sd") sd[j] else mean[j]
      derivative = if (record$type == "sd") c(rep(0, n), ds) else c(dm, rep(0, n))
      system = prior_joint_face_add(system,
        c(value - record$target - record$tolerance, record$target - value - record$tolerance) / unit - face$limit,
        rbind(derivative, -derivative) / unit)
      next
    }
    sign = if (record$operator %in% c("<", "<=")) 1 else -1
    quantiles = c(stats::qnorm(max(0, record$target - record$tolerance)) - face$limit,
      stats::qnorm(min(1, record$target + record$tolerance)) + face$limit)
    if (is.finite(quantiles[2])) system = prior_joint_face_add(system,
      (sign * (record$threshold - mean[j]) - quantiles[2] * sd[j]) / unit,
      matrix(c(-sign * dm, -quantiles[2] * ds) / unit, 1))
    if (is.finite(quantiles[1])) system = prior_joint_face_add(system,
      (quantiles[1] * sd[j] - sign * (record$threshold - mean[j])) / unit,
      matrix(c(sign * dm, quantiles[1] * ds) / unit, 1))
  }
  system
}

# The epigraph objective is unchanged. Probability-bound derivatives now
# include the bilinear compromise * SD term explicitly.
prior_joint_minimax_system = function(x, face) {
  n = length(face$sd)
  p = 2 * n
  local = face
  local$limit = face$coefficient_scale * x[p + 1L]
  system = prior_joint_face_system(head(x, p), local)
  system$value = x[p + 1L]
  system$gradient = c(rep(0, p), 1)
  system$hessian = matrix(0, p + 1L, p + 1L)
  system$jacobian = cbind(system$jacobian, 0)
  for (i in seq_along(system$curvature)) {
    H = matrix(0, p + 1L, p + 1L)
    H[seq_len(p), seq_len(p)] = system$curvature[[i]]
    system$curvature[[i]] = H
  }
  index = system$coefficient_start
  for (record in face$problem$coefficient) {
    probability = record$type == "probability"
    bounds = if (probability) stats::qnorm(c(max(0, record$target - record$tolerance),
      min(1, record$target + record$tolerance))) else c(0, 0)
    rows = index + seq_len(sum(is.finite(bounds)))
    index = index + length(rows)
    derivative = -face$coefficient_scale
    if (probability) {
      j = record$parameter
      derivative = derivative * face$sd[j] * x[n + j] / face$problem$sd[j]
      for (i in rows) {
        system$curvature[[i]][n + j, p + 1L] = -face$coefficient_scale * face$sd[j] / face$problem$sd[j]
        system$curvature[[i]][p + 1L, n + j] = system$curvature[[i]][n + j, p + 1L]
      }
    }
    system$jacobian[rows, p + 1L] = derivative
  }
  prior_joint_face_add(system, -x[p + 1L], matrix(c(rep(0, p), -1), 1))
}

prior_joint_face_evaluate = function(x, face) {
  if (isTRUE(face$minimax)) return(prior_joint_minimax_system(x, face))
  prior_joint_face_system(x, face)
}

# A slack-variable Phase I finds an interior point without enlarging any bands.
prior_joint_face_barrier = function(x, face, mu, phase = FALSE) {
  n = length(face$sd)
  p = 2 * n + as.integer(isTRUE(face$minimax))
  r = x[n + seq_len(n)]
  if (any(!is.finite(x)) || any(r <= 0)) return(list(value = Inf))
  system = prior_joint_face_evaluate(head(x, p), face)
  if (phase) {
    system$value = x[p + 1L]
    system$gradient = c(rep(0, p), 1)
    system$hessian = matrix(0, p + 1L, p + 1L)
    system$g = system$g - x[p + 1L]
    system$jacobian = cbind(system$jacobian, -1)
  }
  if (any(!is.finite(c(system$value, system$g, system$jacobian))) || any(system$g >= 0)) return(list(value = Inf))
  base_gradient = system$gradient
  base_gradient[n + seq_len(n)] = base_gradient[n + seq_len(n)] - mu / r
  gradient = base_gradient - as.vector(crossprod(system$jacobian, mu / system$g))
  base = system$hessian
  for (i in seq_along(system$g)) base[seq_len(p), seq_len(p)] =
    base[seq_len(p), seq_len(p), drop = FALSE] - mu / system$g[i] * system$curvature[[i]]
  diag(base)[n + seq_len(n)] = diag(base)[n + seq_len(n)] + mu / r^2
  weighted = -system$jacobian * sqrt(mu) / system$g
  result = list(value = system$value - mu * sum(log(-system$g)) - mu * sum(log(r)),
    gradient = gradient, hessian = base + crossprod(weighted))
  if (!phase && !isTRUE(face$minimax)) {
    root = tryCatch(chol(base), error = function(error) NULL)
    if (!is.null(root)) result$least_squares = list(matrix = rbind(root, weighted),
      rhs = c(-forwardsolve(t(root), base_gradient), rep(-sqrt(mu), length(system$g))))
  }
  result
}

prior_joint_face_direction = function(system) {
  # Solve the augmented Newton least-squares system directly: forming its
  # normal equations loses the thin face's feasible directions to roundoff.
  if (!is.null(system$least_squares)) {
    A = system$least_squares$matrix
    scale = sqrt(colSums(A^2))
    direction = as.vector(qr.coef(qr(sweep(A, 2, scale, "/"), LAPACK = TRUE), system$least_squares$rhs)) / scale
    if (all(is.finite(direction)) && sum(system$gradient * direction) < 0) return(direction)
  }
  normal = sqrt(pmax(diag(system$hessian), 1e-12))
  spectrum = eigen(system$hessian / outer(normal, normal), symmetric = TRUE)
  direction = as.vector(spectrum$vectors %*%
    (crossprod(spectrum$vectors, -system$gradient / normal) / pmax(spectrum$values, 1e-13))) / normal
  if (any(!is.finite(direction)) || sum(system$gradient * direction) >= 0) {
    return(-system$gradient / pmax(diag(system$hessian), 1))
  }
  direction
}

prior_joint_face_theta = function(x, face) {
  n = length(face$sd)
  result = c(x[seq_len(n)], log(face$sd / face$problem$sd) + log(x[n + seq_len(n)]))
  if (isTRUE(face$minimax)) result = c(result, x[2 * n + 1L])
  result
}

prior_joint_face_search = function(x, face, phase = FALSE) {
  steps = 0L
  p = 2 * length(face$sd) + as.integer(isTRUE(face$minimax))
  for (mu in 10^seq(-2, -11, by = -1)) {
    for (iteration in seq_len(100L)) {
      current = prior_joint_face_barrier(x, face, mu, phase)
      if (!is.finite(current$value)) return(NULL)
      direction = prior_joint_face_direction(current)
      slope = sum(current$gradient * direction)
      rate = 1
      while (rate > 2^-50 && prior_joint_face_barrier(x + rate * direction, face, mu, phase)$value >
          current$value + .01 * rate * slope) rate = rate / 2
      x = x + rate * direction
      steps = steps + 1L
      if (phase) {
        if (max(prior_joint_face_evaluate(head(x, p), face)$g) < 0) {
          return(list(x = head(x, p), iterations = steps))
        }
      } else {
        theta = prior_joint_face_theta(x, face)
        tier = if (isTRUE(face$minimax)) 2L else 3L
        scale = if (isTRUE(face$minimax)) face$coefficient_scale else 1
        original = prior_joint_system(theta, face$problem, tier, face$limit, scale)
        result = list(par = theta, convergence = 0L, iterations = steps,
          objective = original$value, violation = max(c(0, original$g)),
          stationarity = prior_joint_stationarity(original), method = if (isTRUE(face$minimax))
            "Feasible Newton barrier for the minimax epigraph" else "Feasible Newton barrier on the fixed minimax face")
        accurate = !isTRUE(face$minimax) || mu <= 1e-11
        if (accurate && prior_joint_stage_verified(result) && all(prior_joint_details(face$problem, theta)$satisfied)) return(result)
      }
      if (-slope < 1e-13 || rate < 2^-40) break
    }
  }
  NULL
}

# Find a strictly feasible initializer without altering the acceptance bands.
prior_joint_face_interior = function(x, face) {
  if (max(prior_joint_face_evaluate(x, face)$g) < 0) return(x)
  result = prior_joint_face_search(c(x, max(prior_joint_face_evaluate(x, face)$g) + 1), face, TRUE)
  if (is.null(result)) NULL else result$x
}

# Recover a failed KL stage with a more precise minimax value and the same rules.
prior_joint_face_refine = function(start, problem, coefficient_limit, coefficient_scale = 1) {
  if (coefficient_limit <= 0) return(NULL)
  n = length(problem$mean)
  reference_sd = prior_joint_parameters(start, problem)$sd
  face = list(problem = problem, sd = reference_sd,
    limit = coefficient_limit, minimax = TRUE, coefficient_scale = coefficient_scale)
  loss = max(prior_joint_coefficient_metric(start, problem)$violation)
  x = prior_joint_face_interior(c(head(start, n), rep(1, n), loss / coefficient_scale + 1), face)
  if (is.null(x)) return(NULL)
  second = prior_joint_face_search(x, face)
  if (is.null(second)) return(NULL)
  theta = head(second$par, 2 * n)
  limit = max(prior_joint_coefficient_metric(theta, problem)$violation)
  if (limit > coefficient_limit / (1 + 1e-6)) return(NULL)

  # Resolve the existing relative face from a precise minimax value, rather
  # than letting unit-sensitive penalty stopping error set the variance budget.
  # Retain the pre-polishing coordinate scales: a nearly zero minimax SD is
  # a poor KL coordinate scale. This imposes no bound on the selected SDs.
  face = list(problem = problem, sd = reference_sd,
    limit = limit * (1 + 1e-6))
  x = prior_joint_face_interior(c(head(theta, n), rep(1, n)), face)
  if (is.null(x)) return(NULL)
  result = prior_joint_face_search(x, face)
  if (is.null(result) || result$objective > prior_joint_kl(theta, problem)$value +
      1e-8 * max(1, prior_joint_kl(theta, problem)$value)) return(NULL)
  result$minimax = second
  result$coefficient_limit = limit
  result
}
