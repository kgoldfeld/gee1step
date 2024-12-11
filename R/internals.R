# Internal functions

.getD <- function(dd, namesd, rho) {

  # "declare" vars to avoid global NOTE

  v <- NULL
  resid <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  sqrtv <- dd[, sqrt(v)]
  resid <- dd[, resid]
  w <- dd[, w_]
  w1 <- rep(1, times = nrow(d))

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- dv(d, sqrtv, w1)
  sumr <- sum(w * resid)

  return(xrho * sumr * d1)
}

.getW <- function(dd, namesd, rho) {

  # "declare" vars to avoid global NOTE

  v <- NULL

  ###

  d <- as.matrix(dd[, namesd, with = FALSE])
  v <- dd[, v]
  w <- dd[, w_]
  w1 <- rep(1, times = nrow(d))

  sqrtv <- dd[, sqrt(v)]

  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  d1 <- ddv(d, v, w)
  d2_a <- dv(d, sqrtv, w1)
  d2_b <- dv(d, sqrtv, w)

  return( d1 - xrho * d2_a %*% t(d2_b) )
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
  w <- dd[, w_]
  w1 <- rep(1, times = nrow(d))

  # d1 <- dvm(d, resid, w)  # just checking to see if below was a typo in equation (5)
  d1 <- dvm(d, residv, w)
  d2 <- dv(d, sqrtv, w1)
  xrho <- rho / ( (1-rho) + nrow(d) * rho)
  sumr <- sum(w * resid)

  ui <- d1 - xrho * sumr * d2 # from equation (5)

  return(tcrossprod(ui))

}
