#' Predict method for gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param data a required data frame or data.table containing the variables
#' in the model.
#' @param type the type of prediction required. The default ("link") is on the scale of
#' the linear predictors; the alternative "response" is on the scale of the
#' response variable. For this model (a binomial model) the default predictions
#' are of log-odds (probabilities on logit scale) and type = "response" gives
#' the predicted probabilities.
#' @param ... currently disregarded
#' @return a vector of predictions
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' predict(geefit, sampData)
#'
#' @export
predict.gee1step <- function(object, data, type = "link", ...) {

  beta <- object$beta
  X <- object$xnames

  logodds <- as.vector(cbind(1, as.matrix(data[, ..X])) %*% beta)

  ### post-declaration to avoid CHECK note - seems to work
  ..X <- NULL

  if (type == "link") {
    return(logodds)
  } else if (type == "response") {
    return(1/(exp(-logodds) + 1))
  }

}

#' Summarizing gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a \eqn{p \times 4} matrix with columns for the estimated coefficient,
#' its standard error, z-statistic and corresponding (two-sided) p-value.
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' summary(geefit)
#'
#' @export
summary.gee1step <- function(object, ...) {

  se.vb <- sqrt(diag(object$vb))
  z <- object$beta/se.vb
  p.value <- 2*stats::pnorm(-abs(z))

  estimates <- data.frame(est = object$beta, se.err = se.vb, z = z, p.value = p.value)

  rownames(estimates) <-  c("Intercept", labels(stats::terms(object$formula)))

  return(estimates)

}

#' Coef method for gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a vector of length \eqn{p} that represents the parameter estimates.
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' coef(geefit)
#'
#' @export
coef.gee1step <- function(object, ...) {

  beta = as.vector(object$beta)
  names(beta) <- c("(Intercept)", object$xnames)

  return(beta)

}

#' vcov method for gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a matrix of dimension \eqn{p \timts p} that represents the parameter
#' variance-covariance matrix
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' vcov(geefit)
#'
#' @export
vcov.gee1step <- function(object, ...) {

  VCov = object$vb

  names <- c("(Intercept)", object$xnames)
  dimnames(VCov) <- list(names, names)

  return(VCov)
}

#' print method for gee1step model fits
#' @param x a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return nothing is returned
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData, cluster = "site")
#' print(geefit)
#'
#' @export
print.gee1step <- function(x, ...) {

  cat("Call:\n")
  cat(deparse(x$call), "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(coef.gee1step(x))

}



