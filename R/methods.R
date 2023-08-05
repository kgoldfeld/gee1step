#' Summarizing gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a \eqn{p \times 4} matrix with columns for the estimated coefficient,
#' its standard error, z-statistic and corresponding (two-sided) p-value, the
#' ICC, and summary of cluster sizes.
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family="binomial")
#' summary(geefit)
#'
#' @export
summary.gee1step <- function(object, ...) {

  se.vb <- sqrt(diag(object$vb))
  z <- object$beta/se.vb
  p.value <- 2*stats::pnorm(-abs(z))

  estimates <- data.frame(est = object$beta, se.err = se.vb, z = z, p.value = p.value)
  rownames(estimates) <-  object$xnames

  n_clusters <- length(object$cluster_sizes)
  avg_cluster_size <- mean(object$cluster_sizes)
  min_cluster_size <- min(object$cluster_sizes)
  max_cluster_size <- max(object$cluster_sizes)

  clustersum = list(
    n_clusters = n_clusters,
    avg_size = avg_cluster_size,
    min_size = min_cluster_size,
    max_size = max_cluster_size
  )

  sumresult <- list(call = object$call, estimates = estimates, ICC = object$rho, clustersum = clustersum)
  attr(sumresult, "class") <- "summary.gee1step"

  return(sumresult)
}

#' Coef method for gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a vector of length \eqn{p} that represents the parameter estimates.
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family="binomial")
#' coef(geefit)
#'
#' @export
coef.gee1step <- function(object, ...) {

  beta = as.vector(object$beta)
  names(beta) <- object$xnames

  return(beta)

}

#' vcov method for gee1step model fits
#' @param object a fitted model object of class "gee1step".
#' @param ... currently disregarded
#' @return a matrix of dimension \eqn{p \timts p} that represents the parameter
#' variance-covariance matrix
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family="binomial")
#' vcov(geefit)
#'
#' @export
vcov.gee1step <- function(object, ...) {

  VCov = object$vb

  names <- object$xnames
  dimnames(VCov) <- list(names, names)

  return(VCov)
}

#' print method for gee1step model fits
#' @param x a fitted model object of class "gee1step".
#' @param digits the number of significant digits.
#' @param ... currently disregarded
#' @return nothing is returned
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family="binomial")
#' print(geefit)
#'
#'@export
#'
print.gee1step <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(stats::coef(x))) {
    cat("Coefficients:\n")
    print.default(format(stats::coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' print method for summary.gee1step
#' @param x a summary.gee1step object.
#' @param digits the number of significant digits.
#' @param ... currently disregarded
#' @return nothing is returned
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial, cluster = "site", family="binomial")
#' print(geefit)
#'
#' @export
#'
print.summary.gee1step <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  cat("\nCoefficients:\n")
  print(x$estimates)
  cat("\nICC estimate:", x$ICC)
  cat("\nNumber of clusters:", x$clustersum$n_clusters)
  cat("\n\nCluster size summary:\n")
  cat("Avg: ", x$clustersum$avg_size, "  Min: ", x$clustersum$min_size,
      "  Max: ", x$clustersum$max_size, sep = "" )
  cat("\n")
  invisible(x)
}

#' Format to percentage (internal)
#' @param x a number
#' @param ... currently disregarded
#' @return a string with a "%"
#' @noRd
f.perc <- function (x, ...)
  {
  if (length(x)) {
    x <- 100 * x
    ans <- paste0(format(x, trim = TRUE, digits = 3, ...), "%")
    ans
  }
  else character(0)
}

#' Confidence intervals for gee1step model parameters
#' @param object a fitted gee1step model object
#' @param parm a specification of which parameters are to be given confidence
#' intervals, either a vector of numbers or a vector of names. If missing, all
#' parameters are considered.
#' @param level the confidence level required.
#' @param ... currently disregarded
#' @return nothing is returned
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family="binomial")
#' confint(geefit)
#'
#' @export
#'
confint.gee1step <- function(object, parm, level = 0.95, ...) {
  cf <- coef.gee1step(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- f.perc(a)
  fac <- stats::qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov.gee1step(object)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}
