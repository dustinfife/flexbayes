# flexbayes

Minimal R package skeleton for Bayesian extensions to **flexplot**.

## Dependency setup

`DESCRIPTION` lists `flexplot` under `Depends`, so `library(flexbayes)`
attaches `flexplot` and makes its exported functions available:

```r
library(flexbayes)

visualize(model)
estimates(model)
compare.fits(model1, model2)
model.comparison(model1, model2)
```

**Important:** These are **not** exported by `flexbayes` itself; they are
exported by the attached `flexplot` package. Therefore
`flexbayes::estimates` is **not** available. Use `estimates()` after
`library(flexbayes)` or `flexplot::estimates()` explicitly.

## Development

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document("path/to/flexbayes")
devtools::load_all("path/to/flexbayes")
devtools::test("path/to/flexbayes")
```

## Next step

Inspect how each upstream `flexplot` function dispatches. If a function is an
S3 generic, add an appropriate model-specific method and register it in
`NAMESPACE`. Listing `flexplot` in `Depends` does **not** automatically
register Bayesian S3 methods. No Bayesian methods are implemented yet.

Replace the placeholder maintainer email in `DESCRIPTION` before publishing.
