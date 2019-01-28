#' Loglikelihood adjustment of evd fits
#'
#' Loglikelihood adjustment for fitting Generalised Extreme Value (gev) model and
#' Threshold Modelling using generalised Pareto distribution (gpd). The adjustment is
#' based on package \code{evd}.
#'
#' @inherit adj_object params details return references seealso
#'
#' @examples
#' # We need the evd package
#' got_evd <- requireNamespace("evd", quietly = TRUE)
#'
#' if (got_evd) {
#'   library(evd)
#'   # An example from the evd::fgev documentation
#'   uvdata <- evd::rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fgev(uvdata, nsloc = (-49:50)/100)
#'   adj_fgev <- alogLik(M1)
#'   summary(adj_fgev)
#'
#'   # An example from Chandler and Bate (2007)
#'   y <- c(chandwich::owtemps[, "Oxford"], chandwich::owtemps[, "Worthing"])
#'   x <- rep(c(-1, 1), each = length(y) / 2)
#'   owfit <- evd::fgev(y, nsloc = x)
#'   year <- rep(1:(length(y) / 2), 2)
#'   adj_owfit <- alogLik(owfit, cluster = year)
#'   summary(adj_owfit)
#'
#'   # An example from the evd::fpot documentation
#'   uvdata <- evd::rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
#'   M1 <- evd::fpot(uvdata, 1)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#'   # Fit using the pp model, rather than the gpd
#'   M1 <- evd::fpot(uvdata, 1, model = "pp", npp = 365)
#'   adj_fpot <- alogLik(M1)
#'   summary(adj_fpot)
#' }
#' @name evd
NULL
## NULL

#' @rdname evd
#' @export
alogLik.evd <- function(x, cluster = NULL, use_vcov = TRUE, ...) {
  # List of evd objects supported
  supported_by_oolax <- list(gev = c("gev", "uvevd", "evd"),
                             pot = c("pot", "uvevd", "evd"))
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
  if (name_of_class == "gev") {
    if (is.null(x$nsloc)) {
      class(res) <- c("oolax", "chandwich", "evd", "gev", "stat")
    } else {
      class(res) <- c("oolax", "chandwich", "evd", "gev", "nonstat")
    }
  }
  return(res)
}
