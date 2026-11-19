#' Home sales recorded in 2024
#'
#' Sale prices and property characteristics for 294 homes from records
#' retrieved through RentCast's API. Use these data to explore housing-price
#' models or practice expressing beliefs about prices and property effects.
#'
#' @format A data frame with 294 rows and 7 variables:
#' \describe{
#'   \item{sale_date}{Sale date, stored as an R `Date`, from January 2 through
#'   December 31, 2024.}
#'   \item{sale_price}{Recorded sale price in US dollars.}
#'   \item{sqft}{Recorded property area in square feet.}
#'   \item{bedrooms}{Recorded number of bedrooms.}
#'   \item{bathrooms}{Recorded number of bathrooms, including fractional values.}
#'   \item{year_built}{Recorded year of construction.}
#'   \item{age}{Age in years relative to 2024, calculated as
#'   `2024 - year_built`.}
#' }
#' @details
#' Address information (stored under `property_id` in the supplied file), city,
#' and county have been removed. Sale dates have been converted to `Date` values,
#' `square_feet` has been renamed to `sqft`, and `age` has been added.
#' Age uses the fixed reference year 2024, not the current year.
#'
#' Of the 341 supplied records, 47 had a missing value in at least one retained
#' variable and were removed. The 294 remaining records have no missing values.
#' Observed values have not been changed or imputed.
#' @source RentCast API. Records retrieved and supplied by the package author
#'   in `homes_2024.csv`.
#' @examples
#' data(homes_2024)
#' head(homes_2024)
#' summary(homes_2024)
"homes_2024"
