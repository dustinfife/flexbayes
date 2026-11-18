interview_example = function(area_unit=1,currency=1) {
  skip_if_not_installed("rstanarm")
  if (as.character(packageVersion("rstanarm"))!="2.32.2") skip("Verified backend version required")
  d = data.frame(sqft=2014.94+c(-900,-600,-300,300,600,900),age=c(20,40,20,40,30,30),
    bedrooms=c(2,3,4,2,3,4),price=c(210000,230000,250000,290000,310000,330000))
  d$sqft = d$sqft/area_unit; d$price = d$price*currency
  identify_priors(price~sqft+age+bedrooms,d)
}

interview_beliefs = function(area_unit=1,currency=1,second_lower=400000,tolerance=1e-6*currency,probability_tolerance=1e-6) {
  list(update_outcome(list(sqft=1500/area_unit,age=30,bedrooms=3),
      belief(mode=300000*currency,tolerance=tolerance),
      belief(value<250000*currency,probability=.05,tolerance=probability_tolerance)),
    update_outcome(list(sqft=3000/area_unit,age=30,bedrooms=3),
      belief(mode=500000*currency,tolerance=tolerance),
      belief(value<second_lower*currency,probability=.05,tolerance=probability_tolerance)))
}

compatibility_problem = function(p, specifications) {
  id = 0L
  for (i in seq_along(specifications)) {
    specifications[[i]] = unclass(specifications[[i]])
    for (j in seq_along(specifications[[i]]$beliefs)) {
      id = id+1L
      specifications[[i]]$beliefs[[j]]$id = id
    }
  }
  p$outcomes = specifications
  prior_joint_problem(p)
}

test_that("two-house uncertainty ratio conflict is proved using the centered design", {
  p = interview_example(); saved = p
  specs = interview_beliefs(); original = specs
  first = p+specs[[1]]
  expect_warning({bad = first+specs[[2]]},"uncertainty ratio")
  expect_identical(bad$joint$status,"incompatible")
  expect_identical(bad$joint$compatibility$issue,"variance")
  expect_identical(p,saved); expect_identical(specs,original)
  expect_length(bad$outcomes,2)
  expect_identical(bad$outcomes[[1]],first$outcomes[[1]])
  expect_equal(vapply(bad$outcomes[[2]]$beliefs,`[[`,numeric(1),"target"),c(500000,.05))
  expect_identical(bad$selected,p$selected)
  expect_null(bad$joint$stages)
  reference = mean(p$model$data$sqft)
  expect_equal(reference,2014.94)
  expected_X = rbind(c(1,1500-reference,0,0),c(1,3000-reference,0,0))
  expect_equal(unname(bad$joint$design),expected_X)
  cap = (3000-reference)/abs(1500-reference)
  ratio = bad$joint$compatibility$ratio
  expect_equal(ratio$maximum,cap,tolerance=1e-12)
  expect_gt(ratio$required,1.99)
  expect_lt(ratio$maximum,1.92)
  sd = c(50000,100000)/abs(qnorm(.05))
  expect_equal(sd[2]/sd[1],2)
  bounds = bad$joint$compatibility$bounds
  expect_true(all(bounds$Lower_SD<=sd & bounds$Upper_SD>=sd))
  text = paste(capture.output(summary(bad)),collapse="\n")
  expect_match(text,"Conflict concerns: variance",fixed=TRUE)
  expect_match(text,"Configuration 1",fixed=TRUE)
  expect_match(text,"Configuration 2",fixed=TRUE)
  expect_match(text,"beliefs 1, 2",fixed=TRUE)
  expect_match(text,"1. Most likely expected price",fixed=TRUE)
  expect_match(text,"2. P(expected price <",fixed=TRUE)
  expect_match(text,"uncertainty (SD)",fixed=TRUE)
  expect_match(text,"independent Normal",fixed=TRUE)
  expect_match(text,"answers are retained",fixed=TRUE)
  expect_match(text,"Not ready to use",fixed=TRUE)
  expect_false(grepl("N(",text,fixed=TRUE))
  details = paste(capture.output(summary(bad,diagnostics=TRUE)),collapse="\n")
  expect_match(details,"Nonnegative variance ratio bound",fixed=TRUE)
  expect_error(check_prior(bad),"incompatible")
  expect_error(as_rstanarm(bad),"incompatible")
  expect_error(fit(bad),"incompatible")
})

test_that("compatible house uncertainties reconcile and export without approximation", {
  p = interview_example()
  specs = interview_beliefs(second_lower=425000)
  q = p+specs[[1]]+specs[[2]]
  expect_identical(q$joint$status,"compatible")
  expect_null(q$joint$compatibility)
  X = cbind(1,c(1500,3000)-mean(p$model$data$sqft),0,0)
  mu = vapply(q$selected[1:4],`[[`,numeric(1),"mean")
  sd = vapply(q$selected[1:4],`[[`,numeric(1),"sd")
  achieved_mean = as.vector(X%*%mu)
  achieved_sd = sqrt(as.vector(X^2%*%sd^2))
  expect_lte(max(abs(achieved_mean-c(300000,500000))),1e-6)
  expect_lte(max(abs(pnorm(c(250000,425000),achieved_mean,achieved_sd)-.05)),1e-6)
  # Independent two-equation variance solution, not optimizer diagnostics.
  expected = solve(X[,1:2]^2,(c(50000,75000)/abs(qnorm(.05)))^2)
  expect_true(all(expected>0))
  expect_equal(unname(sd[1:2]^2/expected),c(1,1),tolerance=1e-4)
  args = as_rstanarm(q)
  expect_equal(unname(args$prior$location),unname(mu[-1]))
  expect_equal(unname(args$prior$scale),unname(sd[-1]))
  expect_false(args$prior$autoscale)
  set.seed(17)
  check = check_prior(q,100,newdata=data.frame(sqft=c(1500,3000),age=30,bedrooms=3))
  expect_equal(check$expected,check$coefficient_draws%*%t(check$X))
})

test_that("nonnegative variance certificates detect conflicts requiring three configurations", {
  p = joint_example(y~x)
  specs = lapply(seq_along(c(-1,0,2)),function(i) {
    coordinate = c(-1,0,2)[i]+mean(p$model$data$x)
    threshold = 10+qnorm(.05)*sqrt(c(2,1,3)[i])
    update_outcome(list(x=coordinate),belief(mode=10,tolerance=1e-8),
      belief(value<threshold,probability=.05,tolerance=1e-8))
  })
  problem = compatibility_problem(p,specs)
  proof = prior_joint_compatibility(problem)
  expect_identical(proof$issue,"variance")
  expect_match(proof$technical$method,"Separating hyperplane",fixed=TRUE)
  expect_length(proof$configurations,3)
  witness = proof$technical
  expect_true(all(as.vector(crossprod(witness$matrix,witness$weights))>=0))
  positive = witness$weights>0
  supremum = sum(witness$weights[positive]*witness$upper[positive])+
    sum(witness$weights[!positive]*witness$lower[!positive])
  expect_lt(supremum,0)
  # Exact variance relation V(2) = 4 V(-1) - 3 V(0) would require 5,
  # whereas the answers require 3. Every pair separately passes ratio bounds.
  expect_equal(4*2-3*1,5)
  compatible = specs
  compatible[[3]] = update_outcome(list(x=2+mean(p$model$data$x)),belief(mode=10),
    belief(value<10+qnorm(.05)*sqrt(5),probability=.05))
  expect_null(prior_joint_compatibility(compatibility_problem(p,compatible)))
})

test_that("Normal-only probability conflicts and repeated configurations are diagnosed", {
  p = joint_example(y~x)
  statements = list(update_outcome(list(x=3),belief(value<0,probability=.1)),
    update_outcome(list(x=3),belief(value<1,probability=.5)),
    update_outcome(list(x=3),belief(value<4,probability=.9)))
  proof = prior_joint_compatibility(compatibility_problem(p,statements))
  expect_identical(proof$issue,"location and variance")
  expect_setequal(proof$configurations,1:3)
  expect_length(proof$records,3)
  statements[[3]] = update_outcome(list(x=3),belief(value<2,probability=.9))
  expect_null(prior_joint_compatibility(compatibility_problem(p,statements)))
  conflict = list(update_outcome(list(x=3),belief(mode=0)),update_outcome(list(x=3),belief(mean=10)))
  expect_identical(prior_joint_compatibility(compatibility_problem(p,conflict))$issue,"location")
})

test_that("numerical failure remains unresolved and retains every expert response", {
  p = joint_example(); spec = update_outcome(list(x=6,z=1),belief(mean=20))
  local_mocked_bindings(prior_joint_solve=function(...) stop("forced numerical failure"))
  expect_warning({q = p+spec},"[Nn]o mathematical incompatibility")
  expect_identical(q$joint$status,"unresolved")
  expect_null(q$joint$compatibility)
  expect_length(q$outcomes,1)
  expect_equal(q$outcomes[[1]]$beliefs[[1]]$target,20)
  expect_equal(q$outcomes[[1]]$beliefs[[1]]$tolerance,1e-6)
  expect_identical(q$selected,p$selected)
  text = paste(capture.output(summary(q)),collapse="\n")
  expect_match(gsub("[[:space:]]+"," ",text),"does not prove",fixed=TRUE)
  expect_false(grepl("Prior reconciliation: Incompatible",text,fixed=TRUE))
  expect_error(check_prior(q),"unresolved")
  expect_error(as_rstanarm(q),"unresolved")
  expect_error(fit(q),"unresolved")
})

test_that("acceptance intervals and defaults are preserved rather than silently widened", {
  expect_equal(belief(mode=300000)$tolerance,1e-6)
  expect_equal(belief(value<250000,probability=.05)$tolerance,1e-6)
  p = interview_example()
  strict = compatibility_problem(p,interview_beliefs())
  expect_false(is.null(prior_joint_compatibility(strict)))
  approximate = interview_beliefs(tolerance=5000,probability_tolerance=.005)
  problem = compatibility_problem(p,approximate)
  expect_null(prior_joint_compatibility(problem))
  expect_equal(problem$outcome[[1]]$target,300000)
  expect_equal(problem$outcome[[1]]$tolerance,5000)
  expect_equal(problem$outcome[[2]]$tolerance,.005)
  # A generous interval is explicitly different from incompatible point targets.
  expect_equal(prior_joint_moment_bounds(problem)$lower[2] < prior_joint_moment_bounds(strict)$lower[2],TRUE)
})

test_that("variance diagnostics respect equivalent predictor and response units", {
  original = prior_joint_compatibility(compatibility_problem(interview_example(),interview_beliefs()))
  for (units in list(c(1000,1),c(1,.001),c(1000,.001))) {
    changed = prior_joint_compatibility(compatibility_problem(interview_example(units[1],units[2]),
      interview_beliefs(units[1],units[2])))
    expect_identical(changed$issue,original$issue)
    expect_equal(changed$ratio$maximum,original$ratio$maximum,tolerance=1e-10)
    expect_equal(changed$ratio$required,original$ratio$required,tolerance=1e-8)
    expect_equal(changed$bounds$Lower_SD/units[2],original$bounds$Lower_SD,tolerance=1e-8)
  }
})

test_that("known feasible Normal moments are never rejected by necessary-condition checks", {
  p = joint_example(y~x*z+I(x^2))
  at = data.frame(x=c(1,3,6,8),z=c(0,.5,1,0))
  X = model.matrix(delete.response(p$model$terms),at)
  X[,-1] = sweep(X[,-1,drop=FALSE],2,p$model$reference,"-")
  means = as.vector(X%*%c(10,1,-2,.5,.1))
  for (spread in c(.01,1,100)) {
    sd = sqrt(as.vector(X^2%*%(spread*c(1,.3,.2,.1,.05))^2))
    specs = lapply(seq_len(nrow(at)),function(i) {
      threshold = means[i]+qnorm(.05)*sd[i]
      update_outcome(as.list(at[i,]),belief(mean=means[i]),belief(value<threshold,probability=.05))
    })
    expect_null(prior_joint_compatibility(compatibility_problem(p,specs)))
  }
})

test_that("failed compatibility calculations cannot become an impossibility claim", {
  p = joint_example()
  local_mocked_bindings(prior_joint_compatibility=function(...) stop("forced diagnostic failure"),
    prior_joint_solve=function(...) stop("forced solver failure"))
  expect_warning({q = p+update_outcome(list(x=6,z=1),belief(mean=20))},"[Nn]o mathematical incompatibility")
  expect_identical(q$joint$status,"unresolved")
  expect_null(q$joint$compatibility)
  expect_match(q$joint$compatibility_note,"forced diagnostic failure",fixed=TRUE)
  expect_equal(q$outcomes[[1]]$beliefs[[1]]$target,20)
})

test_that("nearby distinct configurations are not declared impossible by numerical rank truncation", {
  p = joint_example(y~x)
  specs = list(update_outcome(list(x=3),belief(mean=0)),
    update_outcome(list(x=3+1e-14),belief(mean=1)))
  problem = compatibility_problem(p,specs)
  expect_false(identical(problem$X[1,],problem$X[2,]))
  # The required slope is huge, but finite: no location contradiction follows.
  expect_null(prior_joint_compatibility(problem))
})
