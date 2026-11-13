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

## Gaussian prior elicitation

The first prior workflow supports full-rank Gaussian identity-link regression
with the verified rstanarm 2.32.2 backend. It does not fit a posterior model.

```r
p = identify_priors(mpg ~ wt + factor(cyl), data = mtcars)
p = update(p, parameter = "wt", belief(mean = -3), belief(sd = 1))
p = update(p, parameter = "wt", remove = 2,
           belief(value < -1, probability = .95))
summary(p)
plot(p) + ggplot2::theme_bw()
```

Defaults are rstanarm's actual defaults, resolved using the observed response
mean and standard deviation and the design-column standard deviations. The
residual standard deviation has the backend's scaled exponential prior.
The centered intercept describes the linear predictor at the **model-matrix
column means** printed by the object. For factors these include contrast-column
means; for interactions they are means of interaction columns, which need not
correspond to a single observed predictor combination. Models without an
intercept are not centered.

Beliefs use exact coefficient names printed by the object, not formula-term
names. Normal and fixed-df Student's t coefficient/intercept candidates are fitted.
Residual-scale elicitation, predictive checks, export,
and posterior fitting are not yet implemented. Unsupported models fail clearly.
No alternative prior is substituted for an unavailable backend default.

Mean, mode, standard deviation, and threshold probabilities are supported.
Targets use an absolute tolerance (default `1e-6`) on their own scale. Acceptable
normal candidates must satisfy every tolerance; their score is the sum of
squared tolerance-standardized residuals. Acceptable candidates rank normal first, then Student's t; see below. Normal slopes can coexist with
the unchanged normal backend defaults.

When beliefs do not uniquely identify a normal distribution, the solution uses
the closest location/scale to the backend default in standardized Euclidean
distance, subject to the elicited linear constraints. This use of the default
is disclosed in the printed object. Invalid positive-scale solutions and
unsatisfied constraints leave beliefs stored with no selected distribution.
This is a normal-family limitation, not proof that no distribution could work.
Displayed belief numbers are resolved before additions; internal IDs remain
stable. Removing all beliefs restores the backend default.

`print(p)` shows a compact table of priors, their sources, and belief counts.
`Elicited*` marks priors that used the backend default to resolve
underdetermination. `summary(p)` adds current belief numbers, warnings, and
centering references; `summary(p, diagnostics = TRUE)` also shows constraint
fit details. Use the displayed per-parameter belief numbers with `remove`.

Probability beliefs are checked for distribution-independent contradictions
before candidate fitting. The check orders `P(value < t)` before
`P(value <= t)` at each threshold, so point masses remain possible. Upper-tail
statements are converted using their exact complements. Equivalent statements
must agree and CDF values must be nondecreasing, allowing each belief's absolute
tolerance. Disjoint regions cannot have combined probability greater than one.
Contradictions retain all beliefs and set the parameter's candidate status to
`contradictory`; summaries identify the conflicting current belief numbers.
Otherwise consistent beliefs without an acceptable normal candidate have status
`no_supported_fit`. Probabilities zero and one are accepted as beliefs, although
proper normal candidates cannot achieve them exactly at finite thresholds.

### Candidate families and compatibility

Coefficient and centered-intercept updates now fit normal and Student's t
candidates. Student's t has **fixed df = 7**; location and scale are fitted,
and its standard deviation is `scale * sqrt(7/5)`. Degrees of freedom are not
optimized or user-configurable in this release. The policy gives finite means
and variances and avoids identifying tail shape from sparse beliefs.
Underdetermined t fits reference the backend default mean and standard
deviation, converting the latter to t scale.

Every candidate must satisfy every belief's absolute tolerance and have positive
scale. Candidates are fitted using the location-scale equations for moments
and quantiles; the resulting original-scale constraints are checked again.
Acceptable candidates rank normal first, then fixed-df t. This explicit
simplicity policy treats all fits within tolerance as acceptable rather than
letting tiny numerical score differences determine tail shape. Scores retain
the sum of squared tolerance-standardized residuals for inspection. Unsupported
or unsuccessful candidates are retained with diagnostics.

`p$candidates$x$alternatives` retains both fits and mathematical ranks;
`p$candidates$x$mathematical` retains the mathematical preference.
`p$selected$x` is the backend-compatible choice. rstanarm requires one shared
family for all slope priors (with vectorized parameters); its centered intercept
has a separate family argument. After every successful update, backend selections are recomputed for the
whole model without refitting unchanged candidates. Backend defaults remain
fixed. Elicited selections are automatic choices among acceptable candidates,
so they may change when another parameter changes. Beliefs, mathematical
preferences, and candidate ranks are preserved. No explicit candidate-selection
API exists yet; a future explicit choice must be treated as a fixed constraint.
The shared slope family maximizes the number of compatible selections, then
minimizes the sum of mathematical ranks, with normal winning remaining ties.
A complete compatible combination is preferred whenever one exists.
If no complete combination exists, blocked parameters remain
`backend_incompatible` with empty selections and overall compatibility false.
Defaults are retained even in this state. Blocked parameters recover
automatically as soon as subsequent revisions allow a compatible combination.
The centered intercept is selected independently. No fitting/export API is implemented yet.

`plot(p, parameter = "x")` compares acceptable candidate densities and labels
which is selected and which is mathematical only. Invalid candidates are not
drawn. `summary(p, diagnostics = TRUE)` reports acceptance and mathematical ranks.
Probability contradictions still skip every candidate fit; otherwise a failure
is reported as `no_supported_fit`, without claiming mathematical impossibility.

Probability thresholds may be variables or numeric expressions evaluated in the
calling environment, for example `belief(value < lower, probability = .01)` or
`belief(value < qt(.01, df = 7), probability = .01)`. Only the right-hand
threshold is evaluated; `value` is a placeholder. The result must be a single
finite real number. Beliefs retain both the resolved numeric threshold for
fitting and the original comparison for display; later changes to variables do
not change an existing belief.
