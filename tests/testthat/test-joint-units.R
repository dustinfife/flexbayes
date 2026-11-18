# First eight Ames records already used by test-joint-outcomes.R. No download.
joint_unit_example = function(area=1,years=1,currency=1,shift=0,offset=0,scaled_tolerance=TRUE) {
  skip_if_not_installed("rstanarm")
  if (as.character(packageVersion("rstanarm"))!="2.32.2") skip("Verified backend version required")
  d = data.frame(sqft=c(1656,896,1329,2110,1629,1604,1338,1280),
    age=c(50,49,52,42,13,12,9,18),bedrooms=c(3,2,3,3,3,3,2,2),
    price=c(215000,105000,172000,244000,189900,195500,213500,191500))
  d$sqft = d$sqft/area+shift
  d$age = d$age/years+shift
  d$price = currency*d$price+offset
  factors = c(currency,currency*area,currency*years,currency)
  tol = if (scaled_tolerance) 1e-6*factors else rep(1e-6,4)
  p = identify_priors(price~sqft+age+bedrooms,d) +
    update_priors("sqft",belief(mode=150*factors[2],tolerance=tol[2]),
      belief(value<100*factors[2],probability=.1)) +
    update_priors("age",belief(mode=-2000*factors[3],tolerance=tol[3]),
      belief(value> -100*factors[3],probability=.01)) +
    update_priors("bedrooms",belief(mode=25000*currency,tolerance=tol[4]),
      belief(value<5000*currency,probability=.1)) +
    update_priors("(Intercept)",belief(mode=250000*currency+offset,tolerance=tol[1]),
      belief(value<1000*currency+offset,probability=.00001))
  at = list(sqft=3000/area+shift,age=30/years+shift,bedrooms=3)
  spec = update_outcome(at,belief(mode=450000*currency+offset,tolerance=1e-6*currency),
    belief(value<350000*currency+offset,probability=.05))
  q = suppressWarnings(p+spec)
  list(p=p,q=q,spec=spec,factors=factors,offset=offset,area=area,years=years,
    shift=shift,currency=currency,at=at)
}

joint_unit_moments = function(fixture, at) {
  p = fixture$q
  weights = c(1,at$sqft/fixture$area+fixture$shift-p$model$reference["sqft"],
    at$age/fixture$years+fixture$shift-p$model$reference["age"],
    at$bedrooms-p$model$reference["bedrooms"])
  mu = vapply(p$selected[1:4],`[[`,numeric(1),"mean")
  sd = vapply(p$selected[1:4],`[[`,numeric(1),"sd")
  c(mean=(sum(weights*mu)-fixture$offset)/fixture$currency,
    sd=sqrt(sum((weights*sd)^2))/fixture$currency)
}

test_that("Ames reconciliation is invariant to area, age and currency units", {
  a = joint_unit_example()
  expect_identical(a$q$joint$status,"reconciled")
  expect_identical(a$q$beliefs,a$p$beliefs)
  expect_identical(a$q$defaults,a$p$defaults)
  expect_identical(a$q$selected$sigma,a$p$selected$sigma)
  expect_true(all(a$q$joint$outcome$satisfied))
  mu = vapply(a$q$selected[1:4],`[[`,numeric(1),"mean")
  sd = vapply(a$q$selected[1:4],`[[`,numeric(1),"sd")
  baseline_sd = vapply(a$p$selected[1:4],`[[`,numeric(1),"sd")
  changes = mu-vapply(a$p$selected[1:4],`[[`,numeric(1),"mean")
  expect_gt(max(abs(changes))/min(abs(changes)),100)
  for (args in list(list(area=1000),list(years=10),list(area=1000,years=10,shift=7),
    list(currency=.001,offset=10),list(area=1000,scaled_tolerance=FALSE))) {
    b = do.call(joint_unit_example,args)
    expect_identical(b$q$joint$status,"reconciled")
    expect_true(all(b$q$joint$outcome$satisfied))
    adjusted = vapply(b$q$selected[1:4],`[[`,numeric(1),"mean")
    adjusted[1] = adjusted[1]-b$offset
    expect_equal(unname((adjusted/b$factors-mu)/baseline_sd),rep(0,4),tolerance=2e-5)
    expect_equal(unname(vapply(b$q$selected[1:4],`[[`,numeric(1),"sd")/b$factors/sd),rep(1,4),tolerance=2e-5)
    expect_identical(b$q$joint$coefficient$satisfied,a$q$joint$coefficient$satisfied)
    expect_equal(b$q$joint$coefficient$normalized_violation,a$q$joint$coefficient$normalized_violation,tolerance=2e-5)
    for (at in list(list(sqft=3000,age=30,bedrooms=3),list(sqft=1500,age=70,bedrooms=2),
      list(sqft=4000,age=5,bedrooms=5))) {
      original = joint_unit_moments(a,at)
      changed = joint_unit_moments(b,at)
      expect_lt(abs(changed[1]-original[1])/original[2],2e-5)
      expect_equal(unname(changed[2]/original[2]),1,tolerance=2e-5)
    }
  }
  moments = joint_unit_moments(a,list(sqft=3000,age=30,bedrooms=3))
  expect_lte(abs(moments[1]-450000),1e-6)
  expect_lte(abs(pnorm(350000,moments[1],moments[2])-.05),1e-6)
})

test_that("dimensionless minimax agrees with a location-only analytic bound", {
  p = joint_example()
  for (name in colnames(p$model$X)) p = p + update_priors(name,belief(mean=0))
  q = suppressWarnings(p+update_outcome(list(x=6,z=1),belief(mean=20)))
  expect_identical(q$joint$status,"reconciled")
  row = c(1,6-mean(p$model$data$x),1-mean(p$model$data$z))
  baseline_sd = vapply(p$selected[1:3],`[[`,numeric(1),"sd")
  # Triangle inequality gives a global lower bound, attained by aligned
  # standardized mean displacements. No optimizer diagnostics enter this bound.
  bound = (20-1e-6-sum(abs(row)*1e-6))/sum(abs(row)*baseline_sd)
  actual = vapply(q$selected[1:3],`[[`,numeric(1),"mean")/baseline_sd
  expect_equal(unname(actual),rep(bound,3),tolerance=2e-5)
  expect_gte(max(q$joint$coefficient$normalized_violation),bound-1e-8)
  expect_lte(max(q$joint$coefficient$normalized_violation),bound*(1+1.01e-6)+1e-9)
  expect_equal(vapply(q$selected[1:3],`[[`,numeric(1),"sd"),baseline_sd,tolerance=1e-10)
  expect_identical(suppressWarnings(p+update_outcome(list(x=6,z=1),belief(mean=20))),q)
})

test_that("probability sacrifices are probit distances beyond actual acceptance bands", {
  p = joint_example(y~1) + update_priors("(Intercept)",belief(mean=0),belief(sd=1),
    belief(value<qnorm(1e-9),probability=1e-9,tolerance=1e-12))
  q = suppressWarnings(p+update_outcome(list(),belief(mode=1)))
  expect_identical(q$joint$status,"reconciled")
  selected = q$selected$`(Intercept)`
  achieved = pnorm(qnorm(1e-9),selected$mean,selected$sd)
  z = (qnorm(1e-9)-selected$mean)/selected$sd
  lower = qnorm(1e-9-1e-12); upper = qnorm(1e-9+1e-12)
  row = q$joint$coefficient[3,]
  expect_equal(row$achieved,achieved)
  expect_equal(row$normalized_discrepancy,z-qnorm(1e-9))
  expect_equal(row$normalized_violation,max(0,lower-z,z-upper))
  expect_identical(row$satisfied,abs(achieved-1e-9)<=1e-12)
  expect_equal(row$discrepancy,achieved-1e-9)
  # If an absolute probability tolerance reaches 0/1, the acceptance band is
  # one-sided on the probit scale, without clipping the target probability.
  baseline = p; baseline$outcomes = list(unclass(update_outcome(list(),belief(mean=1))))
  baseline$beliefs$`(Intercept)`[[3]]$tolerance = 1e-6
  problem = prior_joint_problem(baseline)
  metric = prior_joint_coefficient_metric(c(1,0),problem)
  expect_identical(metric$lower[3],-Inf)
  expect_equal(metric$violation[3],0)
  system = prior_joint_system(c(1,0,1),problem,2L)
  expect_true(all(is.finite(c(system$g,system$jacobian))))
})

test_that("KL preference accounts for both mean and variance changes in the stated direction", {
  fixture = joint_unit_example()
  q = fixture$q; p = fixture$p
  mu = vapply(q$selected[1:4],`[[`,numeric(1),"mean")
  sd = vapply(q$selected[1:4],`[[`,numeric(1),"sd")
  mu0 = vapply(p$selected[1:4],`[[`,numeric(1),"mean")
  sd0 = vapply(p$selected[1:4],`[[`,numeric(1),"sd")
  kl = sum(log(sd0/sd)+(sd^2+(mu-mu0)^2)/(2*sd0^2)-.5)
  reverse = sum(log(sd/sd0)+(sd0^2+(mu-mu0)^2)/(2*sd^2)-.5)
  expect_equal(q$joint$kl,kl,tolerance=1e-12)
  expect_gt(abs(kl-reverse),.01)
  second = head(q$joint$stages$coefficient$par,8)
  expect_lte(kl,prior_joint_kl(second,q$joint$problem)$value+1e-8)
  expect_lte(max(q$joint$coefficient$normalized_violation),q$joint$coefficient_limit*(1+1.01e-6)+1e-9)
  expect_true(all(q$joint$outcome$satisfied))
  # In a convex location-only benchmark, changed SD starts converge to the
  # same analytic KL projection. General probability problems remain local.
  initial = joint_example()+update_outcome(list(x=6,z=1),belief(mean=20))
  problem = initial$joint$problem
  a = prior_joint_solve(rep(0,6),problem,tier=1L)
  b = prior_joint_solve(c(rep(.2,3),rep(.3,3)),problem,tier=1L)
  expect_equal(a$par,b$par,tolerance=1e-6)
  expect_equal(a$objective,b$objective,tolerance=1e-7)
})

test_that("reconciled priors preserve compatible anchors, diagnostics, sampling and export", {
  p = joint_example()+update_priors("x",belief(mean=0))
  at = list(x=mean(p$model$data$x),z=mean(p$model$data$z))
  q = p+update_outcome(at,belief(mean=p$selected$`(Intercept)`$mean))
  expect_equal(q$selected,p$selected,tolerance=1e-12)
  expect_true(all(q$joint$coefficient$satisfied))
  fixture = joint_unit_example(); q = fixture$q
  exported = as_rstanarm(q)
  expect_equal(unname(exported$prior$location),unname(vapply(q$selected[2:4],`[[`,numeric(1),"mean")))
  expect_equal(unname(exported$prior$scale),unname(vapply(q$selected[2:4],`[[`,numeric(1),"sd")))
  expect_identical(q$beliefs,fixture$p$beliefs)
  expect_false(exported$prior$autoscale)
  expect_equal(unname(exported$prior_intercept$location),q$selected$`(Intercept)`$mean)
  set.seed(102)
  check = check_prior(q,10000,newdata=as.data.frame(fixture$at))
  theoretical = joint_unit_moments(fixture,fixture$at)
  expect_lt(abs(mean(check$expected)-theoretical[1])/theoretical[2],.04)
  expect_lt(abs(sd(as.vector(check$expected))/theoretical[2]-1),.04)
  expect_equal(check$expected,check$coefficient_draws%*%t(check$X))
  expect_output(print(summary(q,diagnostics=TRUE)),"baseline-SD.*probit")
  expect_output(print(summary(q,diagnostics=TRUE)),"discrepancy / numerical tolerance")
})

test_that("dimensionless coefficient gradients and epigraph constraints agree with finite differences", {
  fixture = joint_unit_example()
  problem = fixture$q$joint$problem
  theta = fixture$q$joint$theta
  numerical = sapply(seq_along(theta),function(i) {
    step = rep(0,length(theta)); step[i] = 1e-5
    (prior_joint_coefficient_metric(theta+step,problem)$value-
      prior_joint_coefficient_metric(theta-step,problem)$value)/2e-5
  })
  expect_equal(prior_joint_coefficient_metric(theta,problem)$gradient,unname(numerical),tolerance=1e-7)
  theta = c(theta,.4)
  numerical = sapply(seq_along(theta),function(i) {
    step = rep(0,length(theta)); step[i] = 1e-5
    (prior_joint_system(theta+step,problem,2L)$g-
      prior_joint_system(theta-step,problem,2L)$g)/2e-5
  })
  expect_equal(prior_joint_system(theta,problem,2L)$jacobian,unname(numerical),tolerance=1e-7)
})

test_that("standardized beliefs retain entered-unit diagnostic discrepancies", {
  p = joint_example(y~x) + update_priors("x",belief(mean=.5,scale="standardized"),
    belief(sd=.1,scale="standardized")) +
    update_priors("(Intercept)",belief(mean=0),belief(sd=1))
  q = suppressWarnings(p+update_outcome(list(x=6),belief(mean=30)))
  expect_identical(q$joint$status,"reconciled")
  factor = p$beliefs$x[[1]]$standardization_factor
  record = q$joint$coefficient[q$joint$coefficient$id==p$beliefs$x[[1]]$id,]
  expect_equal(record$target,.5)
  expect_equal(record$achieved,q$selected$x$mean/factor)
  expect_equal(record$discrepancy,(q$selected$x$mean-p$beliefs$x[[1]]$target)/factor)
  expect_equal(record$tolerance,1e-6)
  expect_equal(record$normalized_discrepancy,
    record$discrepancy/(p$selected$x$sd/factor))
  expect_identical(q$beliefs,p$beliefs)
})

test_that("outcome conflict certificates do not depend on currency units", {
  for (unit in c(1,.001,1000)) {
    p = joint_example()
    p$model$data$y = p$model$data$y*unit
    p = identify_priors(y~x+z,p$model$data)
    q = p+update_outcome(list(x=6,z=1),belief(mean=20*unit,tolerance=1e-6*unit))
    expect_warning({conflict = q+update_outcome(list(x=6,z=1),
      belief(mean=10*unit,tolerance=1e-6*unit))},"incompatible")
    expect_identical(conflict$joint$status,"incompatible")
    # Multiple beliefs at a configuration must not confuse record indices
    # with configuration indices when scaling the conflict certificate.
    p$outcomes = list(
      unclass(update_outcome(list(x=6,z=1),belief(mean=20*unit),
        belief(value<15*unit,probability=.05))),
      unclass(update_outcome(list(x=1,z=0),belief(mean=0,tolerance=1e-6*unit))),
      unclass(update_outcome(list(x=1,z=0),belief(mean=10*unit,tolerance=1e-6*unit))))
    expect_match(prior_joint_conflict(prior_joint_problem(p)),"incompatible")
  }
})

test_that("zero-sacrifice compatibility and unused backend defaults survive unit changes", {
  for (unit in c(1,1000)) {
    d = data.frame(x=(1:6)/unit+7,z=rep(0:1,3),y=c(2,4,3,7,6,10))
    p = identify_priors(y~x+z,d)+update_priors("x",belief(mode=0,tolerance=1e-6*unit))
    q = p+update_outcome(list(x=6/unit+7,z=.5),belief(mean=20))
    expect_identical(q$joint$status,"compatible")
    expect_true(all(q$joint$coefficient$satisfied))
    expect_lte(abs(q$selected$x$mean)/unit,1e-6)
    expect_equal(q$selected$z,p$selected$z,tolerance=1e-12)
    expect_identical(q$selected$sigma,p$selected$sigma)
  }
})
