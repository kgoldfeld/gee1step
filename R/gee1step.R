#' Estimate parameters using one-step algorithm
#' @useDynLib gee1step, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param family the distribution family: gaussian, binomial, and poisson
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
gee1step <- function(formula, data, cluster, family, weight=NULL, ...) {

  orig.formula <- formula

  # "declare" vars to avoid global NOTE

  cname_ <- NULL
  Y <- NULL
  N <- NULL
  if(!exists("weight")) weight <- NULL

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

  MM <- stats::model.matrix(formula, data = data)

  Y_ <- all.vars(formula)[1]
  X_ <- colnames(MM)

  if (Y_ %in% X_) {
    stop(paste("Outcome variable", Y_, "cannot also be a predictor"))
  }

  namesd <- paste0("d.", X_)

  dx <- data.table::data.table(MM)
  dx[, cname_ := data[, get(cluster)] ]
  dx[, Y := data[, get(Y_)] ]

  formula <- stats::update(formula, Y ~ .)

  N_clusters <- length(unique(dx[, cname_]))

  ### Call proper family

  if (family == "binomial")  {
    result <- gee1step.binomial(dx, formula, X_, Y_, namesd, N_clusters)
  }
  else if (family == "gaussian") {
    result <- gee1step.gaussian(dx, formula, X_, Y_, namesd, N_clusters, weight = weight)
  }
  else if (family == "poisson") {
    result <- gee1step.poisson(dx, formula, X_, Y_, namesd, N_clusters)
  }

  result <- append(
    list(
      call = match.call(),
      formula = orig.formula,
      family = family,
      outcome = Y_,
      xnames = X_,
      model.data = MM,
      cluster_sizes = as.vector(dx[, .N, keyby = cname_][, N])
    ), result)
  attr(result, "class") <- "gee1step"

  return(result)
}

