indirect_row = function(object, parameter) {
  table = prior_report_table(object)
  table[table$Parameter == parameter, , drop = FALSE]
}
indirect_text = function(object, summary = FALSE, diagnostics = FALSE) {
  paste(capture.output(if (summary) print(summary(object, diagnostics = diagnostics)) else print(object)), collapse = "\n")
}

test_that("one-predictor difference beliefs report their actual source for both backends", {
  for (backend in c("rstanarm", "brms")) {
    skip_if_not_installed(backend)
    p = identify_priors(sale_price ~ sqft + age + bedrooms, homes_2024, package = backend)
    q = p + update_priors(difference = list(
      from = list(sqft = 1700, age = 44, bedrooms = 3),
      to = list(sqft = 2000, age = 44, bedrooms = 3)),
      belief(mode = 50000), belief(value > 90000, probability = .1))
    before = q
    row = indirect_row(q, "sqft")
    expect_identical(row$Source, "Elicited via difference")
    expect_identical(row$Beliefs, "2 indirect")
    for (parameter in c("(Intercept)", "age", "bedrooms", "sigma")) {
      row = indirect_row(q, parameter)
      expect_identical(row$Source, if (backend == "brms") "flexbayes default" else "Backend default")
      expect_identical(row$Beliefs, "0")
    }
    expect_equal(prior_indirect_report(q)$total, 2)
    expect_equal(unname(prior_indirect_report(q)$count), c(0, 2, 0, 0, 0))
    for (summary in c(FALSE, TRUE)) {
      text = indirect_text(q, summary)
      expect_match(text, "Elicited via difference", fixed = TRUE)
      expect_match(text, "2 indirect", fixed = TRUE)
      expect_match(text, "do not sum shared rows", fixed = TRUE)
      expect_false(grepl("Adjusted (default)", text, fixed = TRUE))
    }
    expect_identical(q, before)
    expect_length(q$beliefs$sqft, 0)
  }
})

test_that("shared differences count statements without claiming separate identification", {
  p = joint_example()
  q = p + update_priors(difference = list(from = list(x = 1, z = 0), to = list(x = 4, z = 1)),
    belief(mean = 9), belief(value < 12, probability = .9))
  for (parameter in c("x", "z")) {
    expect_identical(indirect_row(q, parameter)$Source, "Reconciled")
    expect_identical(indirect_row(q, parameter)$Beliefs, "2 indirect")
  }
  expect_identical(indirect_row(q, "(Intercept)")$Beliefs, "0")
  metadata = prior_indirect_report(q)
  expect_equal(metadata$total, 2)
  expect_equal(sum(metadata$count), 4)
  expect_match(indirect_text(q), "2 statement(s) total", fixed = TRUE)
  expect_match(indirect_text(q, TRUE), "Reconciled; 2 indirect", fixed = TRUE)
  expect_false(grepl("Elicited via difference", indirect_text(q), fixed = TRUE))
})

test_that("outcome rows and centered references determine participation", {
  p = joint_example()
  q = p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 20))
  expect_true(all(prior_indirect_report(q)$count[c("(Intercept)", "x", "z")] == 1))
  expect_true(all(prior_report_table(q)$Source[1:3] == "Reconciled"))
  expect_identical(indirect_row(q, "sigma")$Beliefs, "0")
  r = p + update_priors(outcome = list(x = mean(p$model$data$x), z = mean(p$model$data$z)), belief(mean = 20))
  expect_identical(indirect_row(r, "(Intercept)")$Source, "Elicited via outcome")
  expect_equal(prior_indirect_report(r)$count[c("x", "z")], c(x = 0, z = 0))
  expect_true(all(prior_report_table(r)$Source[2:4] == "Backend default"))
})

test_that("direct-only reporting and combined direct and indirect counts stay distinct", {
  p = joint_example() + update_priors("x", belief(mean = 0), belief(sd = 1))
  row = indirect_row(p, "x")
  expect_identical(row$Source, "Elicited")
  expect_equal(row$Beliefs, 2)
  expect_false(grepl("Indirect:", indirect_text(p), fixed = TRUE))
  q = suppressWarnings(p + update_priors(difference = list(from = list(x = 1, z = 0), to = list(x = 4, z = 0)), belief(mean = 6)))
  expect_identical(indirect_row(q, "x")$Source, "Reconciled")
  expect_identical(indirect_row(q, "x")$Beliefs, "2 + 1")
  expect_match(indirect_text(q, TRUE), "Reconciled; 2 + 1", fixed = TRUE)
  expect_match(indirect_text(q), "Mixed counts are direct + indirect", fixed = TRUE)
  expect_identical(q$beliefs, p$beliefs)
})

test_that("several outcome/difference specifications use a general reconciliation label", {
  q = joint_example() +
    update_priors(difference = list(from = list(x = 1, z = 0), to = list(x = 4, z = 0)), belief(mean = 6)) +
    update_priors(outcome = list(x = 4, z = 1), belief(mean = 15))
  metadata = prior_indirect_report(q)
  expect_equal(metadata$count, c(`(Intercept)` = 1, x = 2, z = 1, sigma = 0))
  expect_equal(metadata$total, 2)
  expect_identical(indirect_row(q, "x")$Source, "Reconciled")
  expect_identical(indirect_row(q, "x")$Beliefs, "2 indirect")
  text = indirect_text(q, TRUE)
  expect_match(text, "Difference in expected", fixed = TRUE)
  expect_match(text, "At: x = 4, z = 1", fixed = TRUE)
  expect_match(text, "Reconciled", fixed = TRUE)
})

test_that("participation alone does not replace unchanged default/direct source labels", {
  p = joint_example()
  # The already-selected priors satisfy this expected mean. The row still
  # contains all three coefficients, without requiring any prior adjustment.
  mean = sum(prior_joint_design(p, list(list(at = list(x = 6, z = 1)))) *
    vapply(p$selected[1:3], `[[`, numeric(1), "mean"))
  q = p + update_priors(outcome = list(x = 6, z = 1), belief(mean = mean))
  metadata = prior_indirect_report(q)
  expect_true(metadata$verified)
  expect_true(all(metadata$count[1:3] == 1))
  expect_false(any(metadata$changed))
  expect_true(all(prior_report_table(q)$Source == "Backend default"))
  expect_match(indirect_text(q, TRUE), "No coefficient priors changed.", fixed = TRUE)
  expect_match(indirect_text(q, TRUE), "Unchanged participants: (Intercept), x, z", fixed = TRUE)
  expect_match(indirect_text(q), "not change or independent information", fixed = TRUE)
  p = p + update_priors("x", belief(mean = 0))
  q = p + update_priors(outcome = list(x = 6, z = 1), belief(mean = mean))
  expect_identical(indirect_row(q, "x")$Source, "Elicited*")
  expect_identical(indirect_row(q, "x")$Beliefs, "1 + 1")
})

test_that("failed and stale reconciliation never labels attempted priors as elicited", {
  p = joint_example() + update_priors(outcome = list(x = 6, z = 1), belief(mean = 20))
  q = suppressWarnings(p + update_priors(outcome = list(x = 6, z = 1), belief(mean = 0)))
  expect_true(all(grepl("fallback", prior_report_table(q)$Source[1:3])))
  expect_identical(indirect_row(q, "x")$Beliefs, "2 pending")
  expect_false(prior_indirect_report(q)$verified)
  expect_false(grepl("Elicited via|Reconciled", indirect_text(q)))
  expect_match(indirect_text(q), "Pending indirect:", fixed = TRUE)
  expect_match(indirect_text(q, TRUE), "No verified reconciliation to report.", fixed = TRUE)
  expect_identical(q$selected, p$selected)
  stale = p
  stale$selected$x$mean = stale$selected$x$mean + 1
  expect_identical(indirect_row(stale, "x")$Source, "Unverified")
  expect_identical(indirect_row(stale, "x")$Beliefs, "1 pending")
  expect_match(indirect_text(stale, TRUE), "Needs review", fixed = TRUE)
  expect_error(check_prior(q), "incompatible")
})

test_that("R-squared sigma beliefs are not attributed as indirect coefficient constraints", {
  q = joint_example(y ~ x) + update_priors(difference = list(from = list(x = -1), to = list(x = 1)), belief(mean = 4)) +
    update_priors("sigma", belief(value > .6, probability = .9, scale = "r_squared"))
  before = q
  row = indirect_row(q, "sigma")
  expect_identical(row$Source, "Elicited")
  expect_identical(row$Beliefs, "1")
  expect_identical(indirect_row(q, "x")$Beliefs, "1 indirect")
  expect_match(indirect_text(q, TRUE), "Beliefs about model R-squared", fixed = TRUE)
  expect_identical(q$r_squared$status, "verified")
  expect_identical(q, before)
})

test_that("printing and summary auto-printing leave numerical state and RNG unchanged", {
  p = joint_example() + update_priors(outcome = list(x = 6, z = 1), belief(mean = 20))
  before = p
  set.seed(52)
  seed = .Random.seed
  for (diagnostics in c(FALSE, TRUE)) {
    visible = withVisible(summary(p, diagnostics = diagnostics))
    automatic = capture.output(if (visible$visible) print(visible$value))
    explicit = capture.output(print(summary(p, diagnostics = diagnostics)))
    expect_identical(automatic, explicit)
    expect_equal(sum(grepl("Prior reconciliation:", explicit, fixed = TRUE)), 1)
  }
  invisible(capture.output(print(p)))
  expect_identical(p, before)
  expect_identical(.Random.seed, seed)
  expect_identical(p$selected, before$selected)
  expect_identical(p$joint, before$joint)
  expect_identical(p$candidates, before$candidates)
  text = capture.output(print(p))
  expect_lte(max(nchar(text)), 80)
})

test_that("model-matrix cancellation governs transformed terms and repeated specifications", {
  p = joint_example(y ~ x + I(x^2))
  q = p + update_priors(difference = list(from = list(x = -2), to = list(x = 2)),
    belief(mean = 8), belief(value < 12, probability = .9))
  expect_equal(prior_indirect_report(q)$count, c(`(Intercept)` = 0, x = 2, `I(x^2)` = 0, sigma = 0))
  expect_identical(indirect_row(q, "x")$Source, "Elicited via difference")
  expect_identical(indirect_row(q, "I(x^2)")$Source, "Backend default")
  q = joint_example(y ~ x) +
    update_priors(difference = list(from = list(x = 1), to = list(x = 3)), belief(mean = 4)) +
    update_priors(difference = list(from = list(x = 1), to = list(x = 4)), belief(mean = 6))
  expect_identical(indirect_row(q, "x")$Source, "Reconciled")
  expect_identical(indirect_row(q, "x")$Beliefs, "2 indirect")
  expect_equal(prior_indirect_report(q)$total, 2)
})

test_that("small prior changes never hide a compromised coefficient belief", {
  p = joint_example(y ~ 1) + update_priors("(Intercept)", belief(mean = 0, tolerance = 1e-12), belief(sd = 1))
  q = suppressWarnings(p + update_priors(outcome = list(), belief(mean = 1e-9, tolerance = 1e-12)))
  before = q
  expect_true(isTRUE(all.equal(q$selected$`(Intercept)`, p$selected$`(Intercept)`)))
  expect_false(q$joint$coefficient$satisfied[1])
  expect_identical(indirect_row(q, "(Intercept)")$Source, "Reconciled")
  expect_match(indirect_text(q, TRUE), "Mean value", fixed = TRUE)
  expect_identical(q, before)
})

test_that("indirect source labels do not depend on predictor units", {
  p = joint_example(y ~ x)
  q = p + update_priors(difference = list(from = list(x = 1), to = list(x = 4)), belief(mean = 6))
  d = p$model$data
  d$x = d$x * 1e12
  r = identify_priors(y ~ x, d) + update_priors(
    difference = list(from = list(x = 1e12), to = list(x = 4e12)), belief(mean = 6))
  expect_identical(indirect_row(r, "x")$Source, "Elicited via difference")
  expect_identical(indirect_row(r, "x")$Source, indirect_row(q, "x")$Source)
  expect_identical(indirect_row(r, "x")$Beliefs, indirect_row(q, "x")$Beliefs)
  expect_match(indirect_text(r, TRUE), "Elicited via difference", fixed = TRUE)
})
