# Defer deserializing fitted-model environments until the dataset is requested.
# pkgload sources data scripts eagerly, so keep backend loading in the promise.
delayedAssign("marketing_stanreg", {
  local({
    fixture = system.file("extdata", "marketing_stanreg.rda", package = "flexbayes", mustWork = TRUE)
    environment = new.env(parent = baseenv())
    load(fixture, envir = environment)
    get("marketing_stanreg", envir = environment, inherits = FALSE)
  })
})
