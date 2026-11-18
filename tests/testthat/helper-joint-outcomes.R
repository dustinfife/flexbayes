joint_example = function(formula = y ~ x + z) {
  skip_if_not_installed("rstanarm")
  if (as.character(packageVersion("rstanarm")) != "2.32.2") skip("Verified backend version required")
  identify_priors(formula, data.frame(x=1:6,z=rep(0:1,3),y=c(2,4,3,7,6,10)))
}

