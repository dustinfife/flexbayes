test_that("home sale data retain documented fields without location identifiers", {
  target = new.env()
  utils::data("homes_2024", package = "flexbayes", envir = target)
  homes = target$homes_2024

  expect_s3_class(homes, "data.frame")
  expect_identical(dim(homes), c(294L, 7L))
  expect_identical(names(homes), c("sale_date", "sale_price", "sqft",
    "bedrooms", "bathrooms", "year_built", "age"))
  expect_identical(rownames(homes), as.character(seq_len(nrow(homes))))
  expect_s3_class(homes$sale_date, "Date")
  expect_false(anyNA(homes))
  expect_identical(range(homes$sale_date), as.Date(c("2024-01-02", "2024-12-31")))
  expect_identical(homes$age, 2024 - homes$year_built)
  expect_identical(is.na(homes$age), is.na(homes$year_built))
  expect_true("homes_2024" %in% utils::data(package = "flexbayes")$results[, "Item"])
})
