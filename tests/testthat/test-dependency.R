test_that("flexplot exports the four intended interfaces", {
  expected = c("visualize", "model.comparison", "compare.fits", "estimates")
  expect_true(all(expected %in% getNamespaceExports("flexplot")))
})

test_that("flexbayes does not replace flexplot functions", {
  expected = c("visualize", "model.comparison", "compare.fits", "estimates")
  for (name in expected) {
    expect_false(name %in% getNamespaceExports("flexbayes"))
  }
})
