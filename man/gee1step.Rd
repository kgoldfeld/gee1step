% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gee1step.R
\name{gee1step}
\alias{gee1step}
\title{Estimate parameters using one-step algorithm}
\usage{
gee1step(formula, data, cluster, family, ...)
}
\arguments{
\item{formula}{an object of class "formula": a symbolic description of the model to be fitted.}

\item{data}{a required data frame or data.table containing the variables in the model.}

\item{cluster}{the name of the field that identifies the clusters.}

\item{family}{the distribution family: gaussian, binomial, and poisson}

\item{...}{currently disregarded}
}
\value{
a "gee1step" object
}
\description{
Estimate parameters using one-step algorithm
}
\examples{
geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
  cluster = "site", family = "binomial")
geefit

}
\references{
Lipsitz, S., Fitzmaurice, G., Sinha, D., Hevelone, N., Hu, J.,
& Nguyen, L. L. (2017). One-step generalized estimating equations with large
cluster sizes. Journal of Computational and Graphical Statistics, 26(3), 734-737.
}
