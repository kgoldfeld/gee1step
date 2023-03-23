#' Estimate parameters using one-step algorithm
#' @import data.table
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param ... currently disregarded
#' @references Lipsitz, S., Fitzmaurice, G., Sinha, D., Hevelone, N., Hu, J.,
#' & Nguyen, L. L. (2017). One-step generalized estimating equations with large
#' cluster sizes. Journal of Computational and Graphical Statistics, 26(3), 734-737.
#' @return a "gee1step" object
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' geefit
#'
#' @export
gee1step <- function(formula, data, cluster, ...) {

  # "declare" vars to avoid global NOTE

  p <- NULL
  intz <- NULL
  intz1 <- NULL
  resid <- NULL
  wt_ij <- NULL
  rho_ij <- NULL
  sum_r <- NULL
  uss_r <- NULL
  N <- NULL


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

  Y <- all.vars(formula)[1]

  X <- all.vars(formula, unique = FALSE)[-1]
  .X <- paste0(".", X)
  ..X <- paste0("..", X)

  if ( ! all(data[, get(Y)] %in% c(0, 1)) ) {
    stop(paste("Outcome variable", Y, "must be binary!"))
  }

  if (Y %in% X) {
    stop(paste("Outcome variable", Y, "cannot also be a predictor"))
  }

  ###

  dx <- copy(data)
  glmfit <- stats::glm(formula, data = dx, family = stats::binomial(link = "logit"))

  dx[, p := stats::predict.glm(glmfit, type = "response")]
  dx[, intz := sqrt( p * (1-p) )]

  for (i in seq_along(.X)) {
    dx[, .X[i] := get(X[i]) * intz]
  }

  dx[, resid := (get(Y) - p) / sqrt( p * (1 -p) ) ]

  ### Estimate ICC

  drho <- dx[,
             list(.N,
               sum_r = sum(resid),
               uss_r = sum(resid^2)
             ), keyby = get(cluster)]

  drho[, wt_ij := ( N * (N-1) / 2)]
  drho[, rho_ij := sum_r^2 - uss_r ]
  rho <- drho[, (sum(rho_ij)/2) / sum(wt_ij)]

  ### Estimate betas

  getsum <- function(x, .dx) {

    .dd <- .dx[, list(sum(get(x))), keyby = cluster]
    setnames(.dd, "V1", x)
    .dd[]
  }

  dsum <- dx[, .N, keyby = cluster]

  sumlist <- lapply(c("resid", "intz", .X), function(x) getsum(x, dx))
  dz <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = cluster), x = sumlist)

  dsum <- merge(dz, dsum)

  dm <- lapply(c("intz", .X),
               function(x) dsum[, get(x) * sqrt(rho / ( (1 - rho) + N * rho))] )
  M <- do.call(cbind, dm)

  dvd_add <- t(M) %*% M
  dsum[, resid := resid * rho / ( ( 1 - rho ) + N * rho ) ]

  du <- lapply(c("intz", .X), function(x) dsum[, get(x)*resid])
  U <- do.call(cbind, du)
  addU <- apply(U, 2, sum)

  beta <- stats::coef(glmfit)
  bvcov <- stats::vcov(glmfit)

  dvd <- matrix(solve(bvcov) - dvd_add, dim(bvcov)[1], dim(bvcov)[2])

  beta2 <- beta - (solve(dvd) %*% matrix(addU))

  ### Estimate robust standard errors

  dr <- copy(data)
  dvars <- as.matrix(cbind(1, dr[, cbind(mget(X))]))

  lodds <-  dvars %*% beta2

  dr[, p:= 1/(1 + exp(-lodds))]
  dr[, resid := ( get(Y) - p) / sqrt( p * (1 - p) )]
  dr[, intz := sqrt(p * ( 1 - p ))]
  dr[, intz1 := get(Y) - p]

  for (i in seq_along(.X)) {
    dr[, .X[i] := get(X[i]) * intz]
  }

  for (i in seq_along(..X)) {
    dr[, ..X[i] := get(X[i]) * intz1]
  }

  dsum <- dx[, .N, keyby = cluster]
  sumlist <- lapply(c("resid", "intz", "intz1", .X, ..X), function(x) getsum(x, dr))

  dz <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = cluster), x = sumlist)
  dsum <- merge(dsum, dz)

  dv <- lapply(c("intz", .X),
               function(x) dsum[, get(x) * sqrt(rho / ( (1 - rho) + N * rho))] )
  M <- do.call(cbind, dv)

  dvd_add <- t(M) %*% M

  dsum[, resid := resid * rho / ( ( 1 - rho ) + N * rho ) ]

  dr <- lapply(c("intz", .X), function(x) dsum[, get(x)*resid])
  U <- do.call(cbind, dr)
  addU <- apply(U, 2, sum)

  .D <- c("intz", .X)
  .E <- c("intz1", ..X)

  for (i in seq_along(.D)) {
    dsum[, .D[i] := get(.E[i]) - get(.D[i]) * resid ]
  }

  dsum2 <- as.matrix(dsum[, (cbind(mget(.D)))])

  sumcols <- apply(dsum2, 2, sum)

  adjust <- sumcols %*% t(sumcols) / nrow(dsum2)

  uusq <- t(dsum2) %*% dsum2 - adjust
  uusq <- uusq * nrow(dsum2) / ( nrow(dsum2) - 1)

  vb <- solve(dvd) %*% uusq %*% solve(dvd)
  beta2 <- as.vector(beta2)

  n_clusters <- nrow(drho)
  avg_cluster_size <- drho[, mean(N)]
  min_cluster_size <- drho[, min(N)]
  max_cluster_size <- drho[, max(N)]

  rm(.X, ..X)

  result <- list(beta = beta2,
                 vb = vb,
                 rho = rho,
                 clustersum = list (n_clusters = n_clusters,
                                  avg_size = avg_cluster_size, min_size = min_cluster_size, max_size = max_cluster_size),
                 outcome = Y,
                 formula = formula,
                 xnames = X,
                 call = match.call()
  )

  attr(result, "class") <- "gee1step"

  return(result)
}

