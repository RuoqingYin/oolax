# ============================== ismev::gpd.fit ============================= #
#' logLikVec - Adjusted loglikelihood vector for GPD Distribution fitted by
#' ismev::gpd.fit
#'
#' It provides a vector of the contributions to the independence loglikelihood from
#' individual observations for Peaks Over Threshold modelling. Used by
#' \code{logLikVec()}.
#'
#' @param object Object of class "ismev_gpd"。It must come from
#' \code{\link[ismev]{gpd.fit}}.
#'
#' @param pars parameters for the fitted object. Since the object is of class "ismev_gpd",
#' it is not necessary to input the parameters.
#'
#' @param ... Additional optional arguments. At present no optional arguments are used.
#'
#'
#' @inherit adj_object params details return references seealso
#'
#' @export
logLikVec.ismev_gpd <- function(object, pars = NULL, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  # If the parameter estimates have not been provided in pars then extract
  # them from the fitted object
  if (is.null(pars)) {
    pars <- coef(object)
  }
  n_pars <- length(pars)
  #
  if (object$trans & is.null(object$xdat)) {
    stop("Covariate data are needed.  Refit the model using oolax::oogev.fit")
  }
  if (!object$trans) {
    response_data <- object$data
    # If trans = FALSE then there are no covariates and object$data contains
    # the response data
    sigma <- pars[1]
    xi <- pars[2]
  } else {
    # If trans = TRUE then there are covariates
    response_data <- object$xdat
    # The numbers of parameters for sigma and xi
    reg_pars <- sapply(object$model, length)
    npsc <- reg_pars[1] + 1
    npsh <- reg_pars[2] + 1
    # object$sigmat, object$shmat contain design matrices
    # Values of sigma, xi for each observation
    sigma <- object$siglink(object$sigmat %*%
                              (pars[seq(1, length = npsc)]))
    xi <- object$shlink(object$shmat %*%
                          (pars[seq(npsc + 1, length = npsh)]))
  }
  # Calculate the loglikelihood contributions
  if (any(sigma <= 0)) {
    val <- -Inf
  } else {
    val <- evd::dgpd(response_data, loc = object$threshold,
                     scale = sigma, shape = xi, log = TRUE)
  }
  # Return the usual attributes for a "logLik" object
  attr(val, "nobs") <- nobs(object)
  attr(val, "df") <- n_pars
  class(val) <- "logLikVec"
  return(val)
}

#' @export
nobs.ismev_gpd <- function(object, ...) {
  return(length(object$data))
}

#' @export
coef.ismev_gpd <- function(object, ...) {
  val <- object$mle
  names(val) <- ismev_gpd_names(object)
  return(val)
}

#' @export
vcov.ismev_gpd <- function(object, ...) {
  vc <- object$cov
  dimnames(vc) <- list(ismev_gpd_names(object), ismev_gpd_names(object))
  return(vc)
}

#' @export
logLik.ismev_gpd <- function(object, ...) {
  return(logLik(logLikVec(object)))
}

ismev_gpd_names <- function(x) {
  if (x$trans) {
    if (is.null(colnames(x$ydat))) {
      scale_names <- paste0("scale", c("", x$model[[1]]))
      shape_names <- paste0("shape", c("", x$model[[2]]))
    } else {
      cov_names <- colnames(x$ydat)
      scale_names <- paste0("scale", c("", cov_names[x$model[[1]]]))
      shape_names <- paste0("shape", c("", cov_names[x$model[[2]]]))
    }
    val <- c(scale_names, shape_names)
  } else {
    val <- c("scale", "shape")
  }
  return(val)
}

# See ismev_methods.R for nobs, coef, vcov, logLik methods for class "gev.fit"
