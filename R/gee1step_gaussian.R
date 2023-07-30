#' Estimate parameters using one-step algorithm (with Gaussian distribution)
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
gee1step.gaussian <- function(formula, data, cluster, ...) {

  getD <- function(dq) {

    dqi <- cbind(int = 1, dq[, ..X])

    xrho <- rho / ( (1-rho) + nrow(dqi) * rho)
    dd.a <- dv(as.matrix(dqi), dq[, sqrt(v)])
    sumr <- dq[, sum(resid)]

    return(xrho * sumr * dd.a)
  }

  getW <- function(dq) {

    dqi <- cbind(int = 1, dq[, ..X])

    xrho <- rho / ( (1-rho) + nrow(dqi) * rho)
    dd <- ddv(as.matrix(dqi), dq[, v])
    dd.a <- dv(as.matrix(dqi), dq[, sqrt(v)])

    return( dd - xrho * dd.a %*% t(dd.a) )
  }

  getU <- function(dq) {

    dqi <- cbind(int = 1, dq[, ..X])

    d1 <- dvm(as.matrix(dqi), dq[, residv])
    d2 <- dv(as.matrix(dqi), dq[, sqrt(v)])
    xrho <- rho / ( (1-rho) + nrow(dqi) * rho)
    sumr <- dq[, sum(resid)]

    ui <- d1 - xrho * sumr * d2

    return(ui %*% t(ui))

  }

  Y <- all.vars(formula)[1]
  X <- all.vars(formula, unique = FALSE)[-1]

  dx <- data.table::copy(data)
  N_clusters <- length(unique(dx[, get(cluster)]))

  glmfit <- stats::glm(formula, data = dx, family = stats::gaussian) # specific to dist

  dx[, p := stats::predict.glm(glmfit, type = "response")]
  dx[, v := var(resid(glmfit))] # specific to dist
  dx[, resid := (get(Y) - p) / sqrt(v) ]

  ### Estimate ICC

  drho <- dx[,
             list(
               .N,
               sum_r = sum(resid),
               uss_r = sum(resid^2)
             ), keyby = cluster]

  drho[, wt_ij := ( N * (N-1) / 2)]
  drho[, rho_ij := sum_r^2 - uss_r ]
  rho <- drho[, (sum(rho_ij)/2) / sum(wt_ij)]

  ### Estimate beta

  wi <- lapply(1:N_clusters, function(i) getW(dx[site == i]))
  W <- Reduce("+", wi)

  di <- lapply(1:N_clusters, function(i) getD(dx[site == i]))
  D <- Reduce("+", di)

  beta <- stats::coef(glmfit)
  beta2 <- beta - solve(W) %*% D

  ### Robust se

  dr <- data.table::copy(data)
  dvars <- as.matrix(cbind(1, dr[, cbind(mget(X))]))

  dr[, p:= dvars %*% beta2]
  dr[, v:= var(get(Y) - p)]
  dr[, resid := ( get(Y) - p) / sqrt( v )]
  dr[, residv := ( get(Y) - p) / v ]

  wi <- lapply(1:N_clusters, function(i) getW(dr[site == i]))
  W <- Reduce("+", wi)

  ui <- lapply(1:N_clusters, function(i) getU(dr[site == i]))
  U <- Reduce("+", ui)

  ### Get results

  vb <- solve(W) %*% U %*% solve(W)
  beta2 <- as.vector(beta2)

  result <- list(beta = beta2,
                 vb = vb,
                 rho = rho,
                 cluster_sizes = as.vector(drho[, N]),
                 outcome = Y,
                 formula = formula,
                 xnames = X,
                 call = match.call()
  )

  attr(result, "class") <- "gee1step"

  return(result)
}








