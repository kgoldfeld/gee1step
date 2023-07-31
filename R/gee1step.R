#' Estimate parameters using one-step algorithm
#' @useDynLib gee1step, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param family the distribution family: gaussian, binomial
#' @param ... currently disregarded
#' @references Lipsitz, S., Fitzmaurice, G., Sinha, D., Hevelone, N., Hu, J.,
#' & Nguyen, L. L. (2017). One-step generalized estimating equations with large
#' cluster sizes. Journal of Computational and Graphical Statistics, 26(3), 734-737.
#' @return a "gee1step" object
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family = "binomial")
#' geefit
#'
#' @export
gee1step <- function(formula, data, cluster, family, ...) {

  # "declare" vars to avoid global NOTE

  cname_ <- NULL

  ###

  ### Check arguments

  if ( ! inherits(formula, "formula" ) ) {
    stop("The argument `formula` is not properly specified")
  }

  if  (! is.data.frame(data) ) {
    stop("Data must be a data.frame or data.table")
  }

  data <- data.table::as.data.table(data)

  if ( ! (all(all.vars(formula) %in% names(data))) ) {
    stop("Variables in formula not all in data set")
  }

  if ( ! (cluster %in% names(data)) ) {
    stop(paste("Cluster variable", cluster, "is not in data set"))
  }

  Y_ <- all.vars(formula)[1]
  X_ <- all.vars(formula, unique = FALSE)[-1]
  X_ <- c(".xintercept", X_)

  if (Y_ %in% X_) {
    stop(paste("Outcome variable", Y_, "cannot also be a predictor"))
  }

  namesd <- paste0("d.", X_)

  dx <- data.table::copy(data)
  dx[, cname_ := get(cluster)]

  N_clusters <- length(unique(dx[, cname_]))

  ### Call proper family

  if (family == "binomial")  {
    result <- gee1step.binomial(dx, formula, X_, Y_, namesd, N_clusters)
  }
  else if (family == "gaussian") {
    result <- gee1step.gaussian(dx, formula, X_, Y_, namesd, N_clusters)
  }

  return(result)
}

