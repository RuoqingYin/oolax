# ================================ evd::fgev ================================ #
#' logLikVec - Adjusted loglikelihood vector for GEV Distribution fitted by evd::fgev
#'
#' It provides a vector of the contributions to the independence loglikelihood from
#' individual observations for Generalized Extreme Value distribution. Used by
#' \code{alogLik.evd()}.
#'
#' @param object Object of class "gev". It must come from \code{\link[evd]{fgev}}.
#' Also see Details.
#'
#' @param pars parameters for the fitted object. Since the object is of class "gev", it is
#' not necessary to input the parameters.
#'
#' @param ... Additional optional arguments. At present no optional arguments are used.
#'
#' @inherit adj_object params details return references seealso
#'
#' @examples
#' library(evd)
#' uvdata <- rgev(100, loc = 0.13, scale = 1.1, shape = 0.2)
#' M1 <- fgev(uvdata, nsloc = (-49:50)/100)
#' logLik(M1)
#' logLikVec(M1)
#' @export
#'
logLikVec.gev <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # Extract the parameter estimates
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  mu <- pars[1]
  sigma <- pars[n_pars - 1]
  xi <- pars[n_pars]
  if (n_pars > 3) {
    mu_reg <- pars[2:(n_pars - 2)]
    mu <- mu + as.matrix(object$nsloc) %*% mu_reg
  } #chapter 6 (have a look)
  # Calculate the weighted loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- evd::dgev(object$data, loc = as.vector(mu), scale = sigma,
                     shape = xi, log = TRUE)
  }
  attr(val, "nobs") <- object$n
  attr(val, "df") <- n_pars
  class(val) <- "logLik"
  return(val)
}

#' @export
nobs.gev <- function(object) {
  return(object$n)
}

#' @export
coef.gev <- function(object) {
  return(object$estimate)
}

#' @export
vcov.gev <- function(object) {
  return(object$var.cov)
}




