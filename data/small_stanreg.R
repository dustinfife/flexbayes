# Defer deserializing fitted-model environments until the dataset is requested.
# pkgload sources data scripts eagerly, so keep backend loading in the promise.
delayedAssign("small_stanreg", {
  local({
    fixture = system.file("extdata", "small_stanreg.rda", package = "flexbayes", mustWork = TRUE)
    environment = new.env(parent = baseenv())
    load(fixture, envir = environment)
    get("small_stanreg", envir = environment, inherits = FALSE)
  })
})
