# Internal functions

.getD <- function(dd, namesd, rho, weight) {

  # "declare" vars to avoid global NOTE

  v <- NULL
  resid <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  sqrtv <- dd[, sqrt(v)]
  resid <- dd[, resid]

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- dv(d, sqrtv)
  sumr <- sum(weight * resid) # scaled by weights

  return(xrho * sumr * d1)
}

.getW <- function(dd, namesd, rho, weight) {

  # "declare" vars to avoid global NOTE

  v <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  v <- dd[, v]
  sqrtv <- dd[, sqrt(v)]

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- ddv(d, v, weight)
  d2 <- dv(d, sqrtv)
  d3 <- dv2(d, sqrtv, weight)

  return( d1 - xrho * d2 %*% t(d3) )
}

.getU <- function(dd, namesd, rho, weight) {

  # "declare" vars to avoid global NOTE

  v <- NULL
  resid <- NULL
  residv <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  residv <- dd[, residv]
  resid <- dd[, resid]
  sqrtv <- dd[, sqrt(v)]

  d1 <- dvm(d, residv, weight)
  d2 <- dv(d, sqrtv)
  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  sumr <- sum(resid * weight) # scaled by residuals

  ui <- d1 - xrho * sumr * d2

  return(ui %*% t(ui))

}
