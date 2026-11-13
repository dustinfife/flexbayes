prior_example = function() {
  skip_if_not_installed("rstanarm")
  if (as.character(utils::packageVersion("rstanarm")) != "2.32.2") skip("Verified backend version required")
  identify_priors(y ~ x + g, data.frame(y = c(1, 3, 2, 7, 4, 9), x = 1:6, g = factor(rep(c("a", "b"), 3))))
}
