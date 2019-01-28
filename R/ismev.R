#' Loglikelihood adjustment of ismev fits
#'
#' Loglikelihood adjustment for fitting Generalised Extreme Value (gev) model and
#' Threshold Modelling using generalised Pareto distribution (gpd). The adjustment is
#' based on package \code{ismev}.
#'
#' @inherit adj_object params details return references seealso
#'
#' @examples
#' # We need the ismev package
#' got_ismev <- requireNamespace("ismev", quietly = TRUE)
#'
#' if (got_ismev) {
#'   library(ismev)
#'
#'   # An example from Chandler and Bate (2007)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- as.matrix(rep(c(1, -1), each = length(y) / 2))
#'   owfit <- oogev.fit(y, x, mul = 1, sigl = 1, shl = 1, method = "BFGS" )
#'   year <- rep(1:(length(y) / 2), 2)
#'   adj_owfit <- alogLik(owfit, cluster = year, cadjust = FALSE)
#'   summary(adj_owfit)
#'
#'  # An example from the gpd.fit() documentation
#'  data(rain)
#'  fit <- oogpd.fit(rain, 10)
#'  adj_fit <- alogLik(fit)
#'  summary(adj_fit)
#' }
#' @name ismev
NULL
## NULL

#' @rdname ismev
#' @export
alogLik.gev.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evd objects supported
  supported_by_oolax <- list(ismev_gev = c("gev.fit"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_oolax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_oolax)[which(is_supported)]
  class(x) <- name_of_class
  # Call adj_object to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (name_of_class == "ismev_gev") {
    if (!x$trans) {
      class(res) <- c("oolax", "chandwich", "ismev", "gev", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "ismev", "gev", "nonstat")
    }
  }
  return(res)
}


#' @rdname ismev
#' @export
alogLik.gpd.fit <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evd objects supported
  supported_by_oolax <- list(ismev_gpd = c("gpd.fit"))
  # Does x have a supported class?
  is_supported <- NULL
  for (i in 1:length(supported_by_oolax)) {
    is_supported[i] <- identical(class(x), unlist(supported_by_oolax[i],
                                                  use.names = FALSE))
  }
  if (!any(is_supported)) {
    stop(paste("x's class", deparse(class(x)), "is not supported"))
  }
  # Set the class
  name_of_class <- names(supported_by_oolax)[which(is_supported)]
  class(x) <- name_of_class
  # Call adj_object to adjust the loglikelihood
  res <- adj_object(x, cluster = cluster, use_vcov = use_vcov, ...)
  if (name_of_class == "ismev_gpd") {
    if (!x$trans) {
      class(res) <- c("oolax", "chandwich", "ismev", "gpd", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "ismev", "gpd", "nonstat")
    }
  }
  return(res)
}
