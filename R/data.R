#' Simulated data set
#'
#' A simulated dataset containing a binary outcome data clustered at the site level and
#' three covariates.
#'
#' @format A data.table with 9788 rows and 6 variables:
#' \describe{
#'   \item{id}{individual record identifier}
#'   \item{site}{site indentifer}
#'   \item{x1}{continuous covariate}
#'   \item{x2}{continuous covariate}
#'   \item{x3}{continuous covariate}
#'   \item{y}{binary outcome}
#' }
"sampData_binomial"

#' Simulated data set
#'
#' A simulated dataset containing a continuous outcome data clustered at the site level and
#' three covariates.
#'
#' @format A data.table with 10012 rows and 6 variables:
#' \describe{
#'   \item{id}{individual record identifier}
#'   \item{site}{site indentifer}
#'   \item{x1}{continuous covariate}
#'   \item{x2}{continuous covariate}
#'   \item{x3}{continuous covariate}
#'   \item{y}{normally distributed continuous outcome}
#' }
"sampData_gaussian"
