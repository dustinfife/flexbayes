\name{summary.flexbayes_priors}
\alias{summary.flexbayes_priors}
\title{Summarize prior specifications}
\usage{
\method{summary}{flexbayes_priors}(object, diagnostics = FALSE, ...)
}
\arguments{
\item{object}{A prior specification.}
\item{diagnostics}{Include constraint targets, achieved values, standardized residuals, and candidate fit scores.}
\item{...}{Reserved.}
}
\value{A summary.flexbayes_priors object, displayed by its print method.}
\description{Displays a compact prior table, numbered beliefs, relevant warnings, and centered-intercept reference values. Internal belief IDs are not displayed.}
