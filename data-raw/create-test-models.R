# Run from package root AFTER create-test-data.R.
# These fitted-model fixtures avoid recompiling Stan every time tests run.
# They are deliberately small; do not use them to assess convergence.
# Requires brms and rstanarm; a working Stan toolchain may be needed for brms.
stopifnot(requireNamespace("brms", quietly = TRUE))
stopifnot(requireNamespace("rstanarm", quietly = TRUE))
stopifnot(exists("small_bayes"), exists("marketing_bayes"))

small_brms = brms::brm(
  y ~ x + z + a, data = small_bayes,
  chains = 2, iter = 600, warmup = 300, cores = 2,
  seed = 2323, refresh = 0
)
usethis::use_data(small_brms, overwrite = TRUE)

small_brms_interaction = brms::brm(
  y ~ x * a, data = small_bayes,
  chains = 2, iter = 600, warmup = 300, cores = 2,
  seed = 2324, refresh = 0
)
usethis::use_data(small_brms_interaction, overwrite = TRUE)

marketing_brms = brms::brm(
  revenue ~ advertising, data = marketing_bayes,
  chains = 2, iter = 600, warmup = 300, cores = 2,
  seed = 2325, refresh = 0
)
usethis::use_data(marketing_brms, overwrite = TRUE)

small_stanreg = rstanarm::stan_glm(
  y ~ x + z + a, data = small_bayes,
  chains = 2, iter = 600, warmup = 300, cores = 2,
  seed = 2323, refresh = 0
)
usethis::use_data(small_stanreg, overwrite = TRUE)

marketing_stanreg = rstanarm::stan_glm(
  revenue ~ advertising, data = marketing_bayes,
  chains = 2, iter = 600, warmup = 300, cores = 2,
  seed = 2325, refresh = 0
)
usethis::use_data(marketing_stanreg, overwrite = TRUE)
