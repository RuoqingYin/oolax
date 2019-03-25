# ================================ evd::fpot ================================ #
#' logLikVec - loglikelihood vector for GPD
#'
#' It provides a vector of the contributions to the independence loglikelihood from
#' individual observations for Peaks Over Threshold modelling. Used by
#' \code{alogLik.evd()}
#'
#'
#' @param object Object of class "pot"。It must come from \code{\link[evd]{fpot}}.
#'
#' @param pars parameters for the fitted object. Since the object is of class "pot", it is
#' not necessary to input the parameters.
#'
#' @param ... Additional optional arguments. At present no optional arguments are used.
#'
#'
#'
#' @inherit adj_object params details return references seealso
#'
#' @examples
#' library(evd)
#' uvdata <- rgpd(100, loc = 0, scale = 1.1, shape = 0.2)
#' M1 <- fpot(uvdata, 1)
#' logLik(M1)
#' logLikVec(M1)
#' @export
#'
logLikVec.pot <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # Extract the parameter estimates
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  sigma <- pars[1]
  xi <- pars[n_pars]
  #if (n_pars > 3) {
  #  mu_reg <- pars[2:(n_pars - 2)]
  #  mu <- mu + as.matrix(object$nsloc) %*% mu_reg
  #} #chapter 6 (have a look)
  # We ingore the above because no nsloc is allowed in the model
  # i.e. only two parameters in model
  # Calculate the weighted loglikelihood contributions
  if (sigma <= 0) {
    val <- -Inf
  } else {
    val <- evd::dgpd(object$exceedances,
                     loc = object$threshold ,
                     scale = sigma, shape = xi, log = TRUE)
  }
  attr(val, "nobs") <- length(object$exceedances)
  attr(val, "df") <- n_pars
  class(val) <- "logLik"
  return(val)
}

#' @export
nobs.pot <- function(object) {
  return(length(object$exceedances))
}

#' @export
coef.pot <- function(object) {
  return(object$estimate)
}

#' @export
vcov.pot <- function(object) {
  return(object$var.cov)
}




