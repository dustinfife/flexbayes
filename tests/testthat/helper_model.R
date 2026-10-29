get_test_model = function() {
  readRDS(testthat::test_path("fixtures", "brms_mtcars.rds"))
}