# Defer deserializing fitted-model environments until the dataset is requested.
# pkgload sources data scripts eagerly, so keep backend loading in the promise.
delayedAssign("small_brms", {
  local({
    fixture = system.file("extdata", "small_brms.rda", package = "flexbayes", mustWork = TRUE)
    environment = new.env(parent = baseenv())
    load(fixture, envir = environment)
    get("small_brms", envir = environment, inherits = FALSE)
  })
})
