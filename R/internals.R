# Internal functions

.getD <- function(dd, namesd, rho) {

  # "declare" vars to avoid global NOTE

  v <- NULL
  resid <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  sqrtv <- dd[, sqrt(v)]
  resid <- dd[, resid]

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- dv(d, sqrtv)
  sumr <- sum(resid)

  return(xrho * sumr * d1)
}

.getW <- function(dd, namesd, rho) {

  # "declare" vars to avoid global NOTE

  v <- NULL

  ###


  d <- as.matrix(dd[, namesd, with = FALSE])
  v <- dd[, v]
  sqrtv <- dd[, sqrt(v)]

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- ddv(d, v)
  d2 <- dv(d, sqrtv)

  return( d1 - xrho * d2 %*% t(d2) )
}

.getU <- function(dd, namesd, rho) {

  # "declare" vars to avoid global NOTE

  v <- NULL
  resid <- NULL
  residv <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  residv <- dd[, residv]
  resid <- dd[, resid]
  sqrtv <- dd[, sqrt(v)]

  d1 <- dvm(d, residv)
  d2 <- dv(d, sqrtv)
  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  sumr <- sum(resid)

  ui <- d1 - xrho * sumr * d2

  return(ui %*% t(ui))

}
