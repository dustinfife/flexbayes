joint_summary_text = function(object, diagnostics=FALSE) {
  paste(capture.output(print(summary(object,diagnostics=diagnostics))),collapse="\n")
}

test_that("ordinary reconciliation summaries group changed priors and compromised beliefs", {
  p = joint_example(y~1)+update_priors("(Intercept)",belief(mean=0),belief(sd=1))
  q = suppressWarnings(p+update_outcome(list(),belief(mean=10)))
  before = q
  text = joint_summary_text(q)
  expect_match(text,"Prior reconciliation: Successful",fixed=TRUE)
  expect_match(text,"Outcome: y",fixed=TRUE)
  expect_match(text,"At: intercept only",fixed=TRUE)
  expect_match(text,"Outcome beliefs satisfied: Yes",fixed=TRUE)
  expect_match(text,"Requested.*Achieved")
  expect_match(text,"Coefficient priors and beliefs",fixed=TRUE)
  expect_match(text,"Parameter / Belief.*Before.*After")
  expect_match(text,"N(0,1)",fixed=TRUE)
  expect_false(grepl("Compromised coefficient beliefs",text,fixed=TRUE))
  expect_match(text,"  Mean value",fixed=TRUE)
  expect_false(grepl("SD of value",text,fixed=TRUE))
  expect_false(grepl("Ready to use:",text,fixed=TRUE))
  expect_false(grepl("Original belief statements are retained",text,fixed=TRUE))
  expect_false(grepl("Normal(",text,fixed=TRUE))
  for (technical in c("normalized_violation","stationarity","convergence","KL(",
    "residual","Tier-2","probit","scaled_violation")) expect_false(grepl(technical,text,fixed=TRUE))
  advanced = joint_summary_text(q,TRUE)
  for (technical in c("Technical diagnostics","normalized_violation","stationarity","convergence",
    "KL(adjusted || baseline)","Tier-2","probit")) expect_match(advanced,technical,fixed=TRUE)
  expect_identical(q,before)
  expect_identical(q$beliefs,p$beliefs)
})

test_that("outcome summaries show each configuration with readable values and probabilities", {
  p = joint_example()
  q = p+update_outcome(list(x=6,z=1),belief(mode=20),belief(value<15,probability=.05)) +
    update_outcome(list(x=1,z=0),belief(mean=0))
  text = joint_summary_text(q)
  expect_match(text,"Most likely expected y",fixed=TRUE)
  expect_match(text,"P(expected y < 15)",fixed=TRUE)
  expect_match(text,"5%",fixed=TRUE)
  expect_match(text,"At: x = 6, z = 1",fixed=TRUE)
  expect_match(text,"At: x = 1, z = 0",fixed=TRUE)
  expect_equal(sum(grepl("Outcome beliefs",strsplit(text,"\n",fixed=TRUE)[[1]],fixed=TRUE)),3L)
  expect_match(text,"Coefficient priors and beliefs",fixed=TRUE)
  expect_false(grepl("Ready to use:",text,fixed=TRUE))
  expect_equal(prior_joint_summary_number(450000),"450,000")
  expect_equal(prior_joint_summary_number(.00001,probability=TRUE),"0.001%")
  expect_false(prior_joint_summary_number(1e-10,probability=TRUE)=="0%")
  expect_equal(prior_joint_summary_number(Inf),"Not available")
})

test_that("compatible summaries identify unchanged priors and no compromised beliefs", {
  p = joint_example()
  q = p+update_priors("x",belief(mean=0))
  q = q+update_outcome(list(x=mean(p$model$data$x),z=mean(p$model$data$z)),
    belief(mean=p$selected$`(Intercept)`$mean))
  text = joint_summary_text(q)
  expect_match(text,"Successful",fixed=TRUE)
  expect_match(text,"No coefficient priors changed.",fixed=TRUE)
  expect_false(grepl("  Prior",text,fixed=TRUE))
  expect_false(grepl("Ready to use:",text,fixed=TRUE))
})

test_that("incompatible, unresolved and stale priors are not described as usable", {
  p = joint_example()+update_outcome(list(x=6,z=1),belief(mean=20))
  bad = suppressWarnings(p+update_outcome(list(x=6,z=1),belief(mean=0)))
  text = joint_summary_text(bad)
  expect_match(text,"Prior reconciliation: Incompatible",fixed=TRUE)
  expect_match(text,"Outcome beliefs satisfied: No",fixed=TRUE)
  expect_match(text,"Not verified",fixed=TRUE)
  expect_false(grepl("N(",text,fixed=TRUE))
  expect_match(text,"Not ready to use:",fixed=TRUE)
  unresolved = p
  unresolved$joint$status = "unresolved"
  unresolved$joint$message = "convergence=1; stationarity=99"
  text = joint_summary_text(unresolved)
  expect_match(text,"Prior reconciliation: Not resolved",fixed=TRUE)
  expect_match(text,"No verified reconciliation to report",fixed=TRUE)
  expect_match(text,"Not ready to use:",fixed=TRUE)
  expect_false(grepl("stationarity",text,fixed=TRUE))
  expect_match(joint_summary_text(unresolved,TRUE),"stationarity=99",fixed=TRUE)
  stale = p
  stale$selected$x$mean = stale$selected$x$mean+1
  text = joint_summary_text(stale)
  expect_match(text,"Prior reconciliation: Needs review",fixed=TRUE)
  expect_match(text,"Not ready to use:",fixed=TRUE)
})

test_that("joint summary auto-printing and explicit printing display exactly once", {
  p = joint_example()+update_outcome(list(x=6,z=1),belief(mean=20))
  for (diagnostics in c(FALSE,TRUE)) {
    creation = capture.output({ result = withVisible(summary(p,diagnostics=diagnostics)) })
    expect_length(creation,0)
    expect_true(result$visible)
    automatic = capture.output(if (result$visible) print(result$value))
    explicit = capture.output(print(summary(p,diagnostics=diagnostics)))
    expect_equal(automatic,explicit)
    expect_equal(sum(grepl("Prior reconciliation:",explicit,fixed=TRUE)),1)
  }
})

test_that("one grouped comparison omits unchanged priors and satisfied beliefs", {
  p = joint_example()
  for (parameter in c("(Intercept)","x")) p = p+update_priors(parameter,belief(mode=0),belief(sd=1))
  q = suppressWarnings(p+update_outcome(list(x=6,z=.5),belief(mean=20)))
  before = q
  rows = prior_joint_summary_comparison(q)
  expect_equal(rows[[1]],c("(Intercept)","  Prior","  Most likely value","",
    "x","  Prior","  Most likely value"))
  expect_false(any(grepl("SD of value",rows[[1]],fixed=TRUE)))
  expect_false(any(rows[[1]]=="z"))
  expect_equal(rows$Before[rows[[1]]=="  Prior"],c("N(0,1)","N(0,1)"))
  expect_equal(rows$Before[rows[[1]]=="  Most likely value"],c("0","0"))
  achieved = q$joint$coefficient$achieved[!q$joint$coefficient$satisfied]
  expect_equal(rows$After[rows[[1]]=="  Most likely value"],
    vapply(achieved,prior_joint_summary_number,character(1)))
  text = joint_summary_text(q)
  expect_equal(sum(grepl("Coefficient priors and beliefs",strsplit(text,"\n",fixed=TRUE)[[1]],fixed=TRUE)),1)
  expect_match(text,"Other priors unchanged.",fixed=TRUE)
  for (removed in c("Changes to coefficient priors","Compromised coefficient beliefs",
    "Ready to use:","centering reference","backend defaults","Observation-noise","Normal(mean")) {
    expect_false(grepl(removed,text,fixed=TRUE))
  }
  advanced = joint_summary_text(q,TRUE)
  expect_match(advanced,"z : Normal",fixed=TRUE)
  expect_match(advanced,"sd = 1 [satisfied]",fixed=TRUE)
  expect_match(paste(capture.output(print(q)),collapse="\n"),"Normal(",fixed=TRUE)
  expect_identical(q,before)
})

test_that("grouped comparison preserves narrow-column content and small probabilities", {
  old = options(width=80)
  on.exit(options(old),add=TRUE)
  p = joint_example(y~1)+update_priors("(Intercept)",belief(mean=0),belief(sd=1),
    belief(value<qnorm(.00001),probability=.00001))
  q = suppressWarnings(p+update_outcome(list(),belief(mean=10)))
  text = capture.output(print(summary(q)))
  expect_lte(max(nchar(text)),80)
  rows = prior_joint_summary_comparison(q)
  index = grepl("P(value",rows[[1]],fixed=TRUE)
  expect_true(any(index))
  expect_equal(rows$Before[index],"0.001%")
  expect_false(rows$After[index]=="0%")
  long = paste(rep("parameter",12),collapse="_")
  expect_equal(paste(prior_joint_summary_cell(long,20),collapse=""),long)
  lines = capture.output(prior_joint_summary_row(list(long,"N(-1.23456e+15,1.23456e-08)","N(0,1)"),c(24,24,24)))
  expect_lte(max(nchar(lines)),80)
})
