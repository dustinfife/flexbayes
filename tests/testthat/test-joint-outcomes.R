joint_moments = function(p, at) {
  # Independent model-matrix construction and backend centering.
  X = model.matrix(delete.response(terms(p$model$formula)), as.data.frame(at))
  X[,-1] = sweep(X[,-1,drop=FALSE],2,colMeans(p$model$X[,-1,drop=FALSE]),"-")
  mu = vapply(p$selected[colnames(X)], `[[`, numeric(1), "mean")
  sd = vapply(p$selected[colnames(X)], `[[`, numeric(1), "sd")
  c(mean=sum(X*mu), sd=sqrt(sum((X*sd)^2)))
}

test_that("outcome specifications validate and leave their inputs unchanged", {
  p = joint_example()
  saved = p
  spec = update_outcome(list(x=6,z=1),belief(mean=20))
  original = spec
  q = p + spec
  expect_identical(p,saved)
  expect_identical(spec,original)
  expect_identical(q$joint$status,"compatible")
  expect_true(all(q$joint$outcome$satisfied))
  expect_equal(joint_moments(q,list(x=6,z=1))["mean"],c(mean=20),tolerance=1e-6)
  expect_error(update_outcome(list(x=NA),belief(mean=1)),"finite number")
  expect_error(update_outcome(list(x=1,x=2),belief(mean=1)),"uniquely named")
  expect_error(update_outcome(list(x=1),belief(sd=2)),"mean, mode and probability")
  expect_error(p + update_outcome(list(x=1),belief(mean=1)),"exactly these predictors")
  expect_error(1 + spec,"on the left")
})

test_that("location reconciliation agrees with the Mahalanobis benchmark", {
  p = joint_example()
  X = c(1,6-mean(p$model$data$x),1-mean(p$model$data$z))
  mu = vapply(p$selected[1:3],`[[`,numeric(1),"mean")
  sd = vapply(p$selected[1:3],`[[`,numeric(1),"sd")
  target = 20
  expected = mu+sd^2*X/sum(sd^2*X^2)*(target-sum(X*mu))
  expect_equal(prior_joint_projection(matrix(X,nrow=1),target,mu,sd),expected,tolerance=1e-12)
  q = p + update_outcome(list(x=6,z=1),belief(mean=target,tolerance=1e-9))
  actual = vapply(q$selected[1:3],`[[`,numeric(1),"mean")
  expect_equal(actual,expected,tolerance=1e-8)
  expect_equal(vapply(q$selected[1:3],`[[`,numeric(1),"sd"),sd,tolerance=1e-12)
  expect_true(all(abs(actual-mu)>0))
  expect_true(q$joint$underdetermined)
  # At the reference z, it contributes exactly zero and must remain unchanged.
  r = p + update_outcome(list(x=6,z=.5),belief(mode=20))
  expect_identical(r$selected$z,p$selected$z)
  expect_identical(r$selected$sigma,p$selected$sigma)
})

test_that("probability beliefs change scales and agree with independent Normal calculations", {
  p = joint_example()
  q = p + update_outcome(list(x=6,z=1),belief(mode=20),belief(value<15,probability=.05))
  moments = joint_moments(q,list(x=6,z=1))
  expect_equal(unname(moments[1]),20,tolerance=1e-6)
  expect_lte(abs(pnorm(15,moments[1],moments[2])-.05),1e-6)
  expect_true(any(abs(vapply(q$selected[1:3],`[[`,numeric(1),"sd")-
    vapply(p$selected[1:3],`[[`,numeric(1),"sd"))>1e-4))
  expect_identical(p + update_outcome(list(x=6,z=1),belief(mode=20),belief(value<15,probability=.05)),q)
  upper = p + update_outcome(list(x=6,z=1),belief(value>15,probability=.95))
  m = joint_moments(upper,list(x=6,z=1))
  expect_lte(abs(pnorm(15,m[1],m[2],lower.tail=FALSE)-.95),1e-6)
  set.seed(3); before = .Random.seed
  invisible(p + update_outcome(list(x=6,z=1),belief(mode=20)))
  expect_identical(.Random.seed,before)
})

test_that("multiple configurations preserve their common baseline and detect conflicts", {
  p = joint_example()
  q = p + update_outcome(list(x=6,z=1),belief(mean=20)) +
    update_outcome(list(x=1,z=0),belief(mode=0))
  expect_true(all(q$joint$outcome$satisfied))
  expect_identical(q$joint$baseline,p$selected)
  expect_length(q$outcomes,2)
  expect_equal(unname(joint_moments(q,list(x=1,z=0))[1]),0,tolerance=1e-6)
  expect_warning({bad = q + update_outcome(list(x=6,z=1),belief(mean=0))},"incompatible")
  expect_identical(bad$joint$status,"incompatible")
  expect_error(check_prior(bad),"incompatible")
  expect_error(as_rstanarm(bad),"incompatible")
  expect_error(fit(bad),"incompatible")
  expect_output(print(summary(bad,diagnostics=TRUE)),"incompatible")
  expect_warning({prob = p + update_outcome(list(x=6,z=1),
    belief(value<0,probability=.9),belief(value<1,probability=.1))},"probabilit|CDF|monotonic")
  expect_identical(prob$joint$status,"incompatible")
})

test_that("explicit coefficient beliefs take precedence over defaults", {
  p = joint_example() + update_priors("x",belief(mode=0))
  q = p + update_outcome(list(x=6,z=1),belief(mean=20))
  expect_identical(q$joint$status,"compatible")
  expect_lte(abs(q$selected$x$mean),1e-6)
  expect_true(abs(q$selected$`(Intercept)`$mean-p$selected$`(Intercept)`$mean)>1)
  expect_identical(q$beliefs,p$beliefs)
  p = joint_example(y~1) + update_priors("(Intercept)",belief(mean=0),belief(sd=1))
  expect_warning({q = p + update_outcome(list(),belief(mean=10))},"relaxed")
  expect_identical(q$joint$status,"reconciled")
  expect_false(q$joint$coefficient$satisfied[1])
  expect_true(q$joint$coefficient$satisfied[2])
  expect_identical(q$beliefs,p$beliefs)
  expect_output(print(summary(q,diagnostics=TRUE)),"relaxed / unresolved")
  expect_equal(unname(as_rstanarm(q)$prior_intercept$location),10,tolerance=1e-6)
  expect_identical(q$defaults,p$defaults)
  expect_warning({conflict = q + update_outcome(list(),belief(mean=0))},"incompatible")
  expect_output(print(summary(conflict,diagnostics=TRUE)),"incompatible")
})

test_that("centered transformations and exported priors agree with predictive simulation", {
  p = joint_example(y~x*z+I(x^2))
  q = p + update_outcome(list(x=6,z=1),belief(mean=20),belief(value<15,probability=.05))
  expect_equal(unname(joint_moments(q,list(x=6,z=1))[1]),20,tolerance=1e-6)
  a = as_rstanarm(q)
  expect_equal(unname(a$prior$location),unname(vapply(q$selected[2:5],`[[`,numeric(1),"mean")))
  expect_equal(unname(a$prior$scale),unname(vapply(q$selected[2:5],`[[`,numeric(1),"sd")))
  expect_equal(unname(a$prior_intercept$location),q$selected$`(Intercept)`$mean)
  expect_false(a$prior$autoscale)
  set.seed(10)
  check = check_prior(q,n=10000,newdata=data.frame(x=6,z=1))
  m = joint_moments(q,list(x=6,z=1))
  expect_equal(mean(check$expected),unname(m[1]),tolerance=.2)
  expect_equal(sd(as.vector(check$expected)),unname(m[2]),tolerance=.2)
  expect_equal(check$expected,check$coefficient_draws %*% t(check$X))
  expect_true(sd(as.vector(check$simulated))>sd(as.vector(check$expected)))
  expect_s3_class(plot(check,type="outcome"),"ggplot")
  captured = new.env()
  local_mocked_bindings(prior_rstanarm_fit=function(arguments) {captured$args=arguments; "native"})
  expect_identical(fit(q,chains=1,iter=10),"native")
  expect_identical(captured$args$prior,a$prior)
  altered = q; altered$selected$x$mean = altered$selected$x$mean+1
  expect_error(check_prior(altered),"changed")
})

test_that("unsupported families and selected prior families fail clearly", {
  p = joint_example(y~x)
  p = update(p,parameter="x",belief(mean=0),belief(sd=sqrt(7/5)),
    belief(value<1,probability=pt(1,7)))
  unsupported = p + update_outcome(list(x=6),belief(mean=10))
  expect_identical(unsupported$selected, p$selected)
  expect_identical(unsupported$joint$status, "unresolved")
  expect_match(unsupported$joint$message, "Student's t", fixed = TRUE)
  expect_false(unsupported$joint$attempted)
  expect_error(as_rstanarm(unsupported), "unresolved")
  expect_warning(as_rstanarm(unsupported, allow_fallback = TRUE), "fallback")
  p = identify_priors(mpg~wt+factor(cyl),mtcars)
  expect_error(p + update_outcome(list(wt=3,cyl=6),belief(mean=20)),"numeric predictors")
  p = prior_example()
  expect_error(p + update_outcome(list(x=6,g=1),belief(mean=10)),"numeric predictors")
  p = identify_priors(y~x,data.frame(x=1:6,y=rep(0:1,3)),family=binomial())
  expect_error(p + update_outcome(list(x=6),belief(mean=.5)),"Gaussian identity")
  p = joint_example(y~0+x)
  expect_error(p + update_outcome(list(x=6),belief(mean=10)),"intercept")
  p = joint_example()
  boundary = p + update_outcome(list(x=6,z=1),belief(value<10,probability=0))
  expect_true(boundary$joint$status %in% c("compatible","reconciled"))
  expect_equal(boundary$outcomes[[1]]$beliefs[[1]]$target,0)
  expect_lte(boundary$joint$outcome$achieved,1e-6)
})

test_that("the Ames worked example meets location and probability tolerances in dollars", {
  # First eight records of https://datasci112.stanford.edu/data/housing.tsv
  # Ames housing; age is sale year minus construction year. No network in tests.
  ames = data.frame(sqft=c(1656,896,1329,2110,1629,1604,1338,1280),
    age=c(50,49,52,42,13,12,9,18),bedrooms=c(3,2,3,3,3,3,2,2),
    price=c(215000,105000,172000,244000,189900,195500,213500,191500))
  p = identify_priors(price~sqft+age+bedrooms,ames)
  q = p + update_outcome(list(sqft=3000,age=30,bedrooms=3),
    belief(mode=450000),belief(value<350000,probability=.05))
  m = joint_moments(q,list(sqft=3000,age=30,bedrooms=3))
  expect_lte(abs(m[1]-450000),1e-6)
  expect_lte(abs(pnorm(350000,m[1],m[2])-.05),1e-6)
  expect_true(all(q$joint$outcome$satisfied))
  expect_identical(q$selected$sigma,p$selected$sigma)
})


test_that("reconciled expected curves use the existing model renderer", {
  skip_if_not("prepare_flexplot" %in% getNamespaceExports("flexplot"), "New public flexplot renderer required")
  q = joint_example() + update_outcome(list(x=6,z=1),belief(mean=20))
  set.seed(12)
  check = check_prior(q,n=50)
  p = plot(check,type="model",display="draws",draws=3)
  expect_s3_class(p,"ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})


test_that("overlapping tolerances do not require identical outcome targets", {
  p = joint_example()
  q = p + update_outcome(list(x=6,z=1),belief(mean=20,tolerance=1)) +
    update_outcome(list(x=6,z=1),belief(mode=21.9999,tolerance=1))
  expect_identical(q$joint$status,"compatible")
  expect_true(all(q$joint$outcome$satisfied))
})

test_that("analytic gradients agree with independent finite differences", {
  q = joint_example() + update_outcome(list(x=6,z=1),belief(mode=20),
    belief(value<15,probability=.05))
  theta = q$joint$theta
  numerical = sapply(seq_along(theta),function(i) {
    step = rep(0,length(theta)); step[i] = 1e-5
    (prior_joint_values(theta+step,q$joint$problem)$value-
      prior_joint_values(theta-step,q$joint$problem)$value)/2e-5
  })
  expect_equal(prior_joint_values(theta,q$joint$problem)$gradient,numerical,tolerance=1e-7)
})


test_that("variance allocation minimizes KL rather than dividing uncertainty equally", {
  q = joint_example() + update_outcome(list(x=6,z=1),belief(mean=20),
    belief(value<15,probability=.05))
  problem = q$joint$problem
  achieved = joint_moments(q,list(x=6,z=1))
  a = (problem$X[1,]*problem$sd)^2
  # Independent Lagrange-multiplier solution for fixed achieved variance.
  lambda = uniroot(function(l) sum(a/(1+l*a))-achieved[2]^2,
    c(0,100),tol=1e-12)$root
  expected = problem$sd^2/(1+lambda*a)
  actual = vapply(q$selected[1:3],function(p)p$sd^2,numeric(1))
  expect_equal(actual,expected,tolerance=1e-4)
  expect_true(all(vapply(q$joint$stages,function(s)s$stationarity<1e-4,logical(1))))
  expect_true(all(vapply(q$joint$stages,function(s)s$violation<1e-8,logical(1))))
})

test_that("coefficient updates rebuild the anchor and retain all outcome specifications", {
  p = joint_example()
  q = p + update_outcome(list(x=6,z=1),belief(mean=20))
  r = q + update_priors("z",belief(mean=2))
  independent = p + update_priors("z",belief(mean=2))
  expect_identical(r$joint$baseline,independent$selected)
  expect_identical(r$outcomes,q$outcomes)
  expect_true(all(r$joint$outcome$satisfied))
  expect_true(all(r$joint$coefficient$satisfied))
})
