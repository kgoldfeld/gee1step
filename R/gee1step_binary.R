# Estimate parameters using one-step algorithm (with Gaussian distribution)
# @param formula an object of class "formula": a symbolic description of the model to be fitted.
# @param data a required data frame or data.table containing the variables in the model.
# @param cluster the name of the field that identifies the clusters.
# @param ... currently disregarded
# @references Lipsitz, S., Fitzmaurice, G., Sinha, D., Hevelone, N., Hu, J.,
# & Nguyen, L. L. (2017). One-step generalized estimating equations with large
# cluster sizes. Journal of Computational and Graphical Statistics, 26(3), 734-737.
# @return a "gee1step" object
# @examples
# geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
# geefit
#
gee1step.binomial <- function(dx, formula, X, Y, namesd, cluster, N_clusters,...) {

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

  ###

  dx[, .xintercept := 1]
  dr <- data.table::copy(dx) # for robust se

  glmfit <- stats::glm(formula, data = dx, family = stats::binomial) # specific to dist

  dx[, p := stats::predict.glm(glmfit, type = "response")]
  dx[, v := p*(1-p)] # specific to dist
  dx[, resid := (get(Y) - p) / sqrt(v) ]

  dX <- dx[, ..X]
  dX <- dX * dx[, p*(1-p)] # specific to dist
  setnames(dX, namesd)
  dx <- cbind(dx, dX)

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

  wi <- lapply(1:N_clusters, function(i) .getW(dx[cname_ == i], namesd, rho))
  W <- Reduce("+", wi)

  di <- lapply(1:N_clusters, function(i) .getD(dx[cname_ == i], namesd, rho))
  D <- Reduce("+", di)

  beta <- stats::coef(glmfit)
  beta2 <- beta - solve(W) %*% D

  ### Robust se

  dvars <- as.matrix(dr[, ..X])

  dr[, lodds := dvars %*% beta2] # specific to dist
  dr[, p := 1/(1 + exp(-lodds))] # specific to dist
  dr[, v := p * (1 - p)] #  # specific to dist
  dr[, resid := ( get(Y) - p) / sqrt( v )]
  dr[, residv := ( get(Y) - p) / v ]

  dR <- dr[, ..X]
  dR <- dR * dr[, p*(1-p)] # specific to dist
  setnames(dR, namesd)
  dr <- cbind(dr, dR)

  wi <- lapply(1:N_clusters, function(i) .getW(dr[cname_ == i], namesd, rho))
  W <- Reduce("+", wi)

  ui <- lapply(1:N_clusters, function(i) .getU(dr[cname_ == i], namesd, rho))
  U <- Reduce("+", ui)

  ### Get results

  vb <- solve(W) %*% U %*% solve(W)


  result <- list(beta =  as.vector(beta2),
                 vb = vb,
                 rho = rho,
                 cluster_sizes = as.vector(drho[, N]),
                 outcome = Y,
                 formula = formula,
                 xnames = X[-1],
                 family = "binomial",
                 call = match.call()
  )

  attr(result, "class") <- "gee1step"

  return(result)
}








