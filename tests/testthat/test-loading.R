test_that("loading does not initialize optional Bayesian backends", {
  skip_if_not_installed("pkgload")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstanarm")
  path = normalizePath(system.file(package = "flexbayes"))
  script = tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    paste0("path = ", encodeString(path, quote = '"')),
    'library(flexplot)',
    'if (file.exists(file.path(path, "Meta", "package.rds"))) {',
    '  library(flexbayes, lib.loc = dirname(path))',
    '} else {',
    '  pkgload::load_all(path, helpers = FALSE, attach_testthat = FALSE)',
    '}',
    'stopifnot(!any(c("brms", "rstanarm", "rstan", "StanHeaders") %in% loadedNamespaces()))',
    'for (class in c("stanreg", "brmsfit")) {',
    '  stopifnot(is.function(getS3method("estimates", class)))',
    '  stopifnot(is.function(getS3method("visualize", class)))',
    '}',
    'stopifnot(inherits(flexbayes::belief(value < -1, probability = .1), "prior_belief"))',
    'data("small_bayes", package = "flexbayes")',
    'stopifnot(is.data.frame(small_bayes))',
    'stopifnot(!any(c("brms", "rstanarm") %in% loadedNamespaces()))',
    'models = c("marketing_brms", "marketing_stanreg", "small_brms", "small_brms_interaction", "small_stanreg")',
    'for (name in models) {',
    '  expected = if (grepl("stanreg", name)) "stanreg" else "brmsfit"',
    '  stopifnot(inherits(get(name, envir = as.environment("package:flexbayes")), expected))',
    '  target = new.env()',
    '  data(list = name, package = "flexbayes", envir = target)',
    '  stopifnot(inherits(get(name, envir = target), expected))',
    '}',
    'target = new.env()',
    'data(marketing_brms, envir = target)',
    'stopifnot(inherits(target$marketing_brms, "brmsfit"))',
    'discovered = data(package = "flexbayes")$results[, "Item"]',
    'stopifnot(all(models %in% discovered))',
    'cat("LAZY_LOAD_OK\n")'
  ), script)
  output = system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  expect_true(any(grepl("LAZY_LOAD_OK", output)))
})

test_that("fitted datasets remain available on explicit request", {
  skip_if_not_installed("brms")
  skip_if_not_installed("rstanarm")
  environment = new.env()
  for (name in c("small_brms", "small_brms_interaction", "marketing_brms", "small_stanreg", "marketing_stanreg")) {
    data(list = name, package = "flexbayes", envir = environment)
    expected = if (grepl("stanreg", name)) "stanreg" else "brmsfit"
    expect_s3_class(get(name, envir = environment), expected)
  }
})
