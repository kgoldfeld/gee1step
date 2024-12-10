# Internal estimation function
#
gee1step.dist <- function(orig.data, dx, formula, family, X_, Y_, namesd, N_clusters, ...) {

  # "declare" vars to avoid global NOTE

  p <- NULL
  resid <- NULL
  wt_ij <- NULL
  rho_ij <- NULL
  sum_r <- NULL
  uss_r <- NULL
  N <- NULL
  cname_ <- NULL
  v <- NULL
  residv <- NULL
  Y <- NULL

  ###

  dr <- data.table::copy(dx) # for robust se

  xnames <- names(dx)
  xnames <- xnames[1:(length(xnames) - 3)] # exclude w, Y and cluster

  ### Verify that weights should or should not be included in initial model fit:

  if (family == "binomial") {
    bform <- stats::update(formula, cbind(Y, 1-Y) ~ . )
    glmfit <- stats::glm(bform, data = orig.data, family = stats::binomial, weights = w_)
  } else if (family == "poisson") {
    glmfit <- stats::glm(formula, data = orig.data, family = stats::poisson, weights = w_)
  } else if (family == "gaussian") {
    glmfit <- stats::glm(formula, data = orig.data, family = stats::gaussian, weights = w_)
  }

  dx[, p := stats::predict.glm(glmfit, type = "response")]

  if (family == "binomial") {
    dx[, v := p*(1-p)]
  } else if (family == "poisson") {
    dx[, v := p]
  } else if (family == "gaussian") {
    dx[, v := stats::var(resid(glmfit))]
  }


  dx[, resid := (Y - p) / sqrt(v) ]

  dX <- dx[, X_, with = FALSE] # no modification (below) for gaussian

  if (family == "binomial") {
    dX <- dX * dx[, p*(1-p)]
  } else if (family == "poisson") {
    dX <- dX * dx[, p]
  }

  setnames(dX, namesd)
  dx <- cbind(dx, dX)

  ### Estimate ICC -

  drho <- dx[,
             list(
               .N,
               sum_r = sum(resid),
               uss_r = crossprod(resid),
               wgt_m = mean(w_)
             ), keyby = cname_]

  # NEW WEIGHTED

  # drho[, wt_ij := wgt_m * ( N * (N-1) / 2)]
  # drho[, rho_ij := wgt_m * (sum_r^2 - uss_r) ] # added parentheses 12/2/24

  # OLD NO WEIGHTS

  drho[, wt_ij := ( N * (N-1) / 2)]
  drho[, rho_ij := sum_r^2 - uss_r ]

  rho <- drho[, (sum(rho_ij)/2) / sum(wt_ij)]


  # change from 1:N_clusters to unique(dx[, cname_]) in lapply 2/12/24

  wi <- lapply(unique(dx[, cname_]), function(i) .getW(dx[cname_ == i], namesd, rho))
  W <- Reduce("+", wi)

  di <- lapply(unique(dx[, cname_]), function(i) .getD(dx[cname_ == i], namesd, rho))
  D <- Reduce("+", di)

  beta <- stats::coef(glmfit)
  beta2 <- beta - solve(W) %*% D  # from equation (12) in Lipsitz

  ### Robust se

  dvars <- as.matrix(dr[, X_, with = FALSE])
  dr[, p := dvars %*% beta2]

  if (family == "binomial") {
    dr[, p := 1/(1 + exp(-p))]
    dr[, v := p * (1 - p)]
  } else if (family == "poisson") {
    dr[, p := exp(p)]
    dr[, v := p]
  } else if (family == "gaussian") {
    dr[, v:= stats::var( Y - p )]
  }

  dr[, resid := ( Y - p) / sqrt( v )]
  dr[, residv := ( Y - p) / v ]

  dR <- dr[, X_, with = FALSE] # no modification (below) for gaussiaan

  if (family == "binomial") {
    dR <- dR * dr[, p*(1-p)]
  } else if (family == "poisson") {
    dR <- dR * dr[, p]
  }

  setnames(dR, namesd)
  dr <- cbind(dr, dR)

  wi <- lapply(unique(dr[, cname_]), function(i) .getW(dr[cname_ == i], namesd, rho))
  W <- Reduce("+", wi)

  # print(solve(W)) # <------------ delete

  ui <- lapply(unique(dr[, cname_]), function(i) .getU(dr[cname_ == i], namesd, rho))
  U <- Reduce("+", ui)

  # print(U) # <------------ delete

  vb <- solve(W) %*% U %*% solve(W)

  # print(vb) # <------------ delete

  ### Get results

  result <- list(beta = as.vector(beta2),
                 vb = vb,
                 rho = rho
  )

  return(result)
}

#' Estimate parameters using one-step algorithm
#' @useDynLib gee1step, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param family the distribution family: gaussian, binomial, and poisson
#' @param weights optional name of the field that identifies observation weights to be applied in the
#' regression, allowing differential importance to be assigned to individual observations.
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
gee1step <- function(formula, data, cluster, family, weights = NULL, ...) {

  orig.formula <- formula

  # "declare" vars to avoid global NOTE

  cname_ <- NULL
  Y <- NULL
  N <- NULL

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

  if ( !is.null(weights)) {
    if (!(weights%in% names(data)) ) {
      stop(paste("Weights variable", weights, "is not in data set"))
    }
  }

  MM <- stats::model.matrix(formula, data = data)

  Y_ <- all.vars(formula)[1]
  X_ <- colnames(MM)

  orig.data <- copy(data) # need original data set to fit glm
  setnames(orig.data, Y_, "Y") # since formula outcome changes (below)
  if (is.null(weights)) {
    orig.data[, w_ := 1]
  } else {
    orig.data[, w_ := data[, get(weights)]]
  }

  if (Y_ %in% X_) {
    stop(paste("Outcome variable", Y_, "cannot also be a predictor"))
  }

  namesd <- paste0("d.", X_)

  dx <- data.table::data.table(MM)
  if (is.null(weights)) {
    dx[, w_ := 1]
  } else dx[, w_ := data[, get(weights)] ]

  dx[, cname_ := data[, get(cluster)] ]
  dx[, Y := data[, get(Y_)] ] # changing outcome to generic "Y"
  formula <- stats::update(formula, Y ~ .)

  N_clusters <- length(unique(dx[, cname_]))

  mod.fit <- gee1step.dist(orig.data, dx, formula, family, X_, Y_, namesd, N_clusters)

  result <- append(
    list(
      call = match.call(),
      formula = orig.formula,
      family = family,
      outcome = Y_,
      xnames = X_,
      model.data = MM,
      weights= weights,
      cluster_sizes = as.vector(dx[, .N, keyby = cname_][, N])
    ), mod.fit)
  attr(result, "class") <- "gee1step"

  return(result)
}

