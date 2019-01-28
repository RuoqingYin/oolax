#' Loglikelihood adjustment of fitted model objects
#'
#' Loglikelihood adjustment of fitted model objects
#'
#' @param x A fitted model object.
#'
#' @param cluster see details in \code{\link[chandwich]{adjust_loglik}}.
#'
#' @param use_vcov A logical scalar. If \code{use_vcov = TRUE} and the method for
#'        \code{x} exists, we use it to estimate the Hessian of the independence
#'        loglikelihood and pass it as \code{H} to \code{\link[chandwich]{adjust_loglik}}.
#'        Otherwise, set \code{H} as \code{NULL}.
#'
#' @param use_sandwich A logical scalar. If \code{use_sandwich = TRUE}, we use
#'        \code{\link[sandwich]{meat}} or \code{\link[sandwich]{meatCL}} to estimate the
#'        \code{V} in \code{\link[chandwich]{adjust_loglik}}. Otherwise, set \code{V} as
#'        \code{NULL}.
#'
#' @param ... Further arguments to be passed to the functions in the
#'   sandwich package \code{\link[sandwich]{meat}}, if \code{cluster = NULL},
#'   or \code{\link[sandwich:vcovCL]{meatCL}}, otherwise.
#'
#'
#' @return An object of class \code{c("oolax", "chandwich")} with the same structure as an
#'         object returned from \code{\link[chandwich]{adjust_loglik}}.
#'
#'
#'
#' @export
adj_object <- function(x, cluster = NULL, use_sandwich = TRUE,
                       use_vcov = TRUE, ...) {
  # Find all available methods for x
  find_methods_fn <- function(i) as.vector(utils::methods(class = class(x)[i]))
  all_methods <- unlist(sapply(1:length(class(x)), find_methods_fn))
  #
  # Set logLikVec (if a method exists) ----------
  #
  # NB: %in% is value matching and "returns a vector of the positions of (first)
  #     matches of its first argument in its second
  has_logLikVec_method <- paste0("logLikVec.", class(x)) %in% all_methods
  if (!any(has_logLikVec_method)) {
    stop("A logLikVec method must be available for x")
  }
  loglik_fn <- function(pars, fitted_object) {
    return(logLikVec(fitted_object, pars = pars))
  }
  #
  # Set H, but not if use_cov = FALSE or no vcov method exists ----------
  #
  if (!use_vcov) {
    H <- NULL
  } else {
    # Check whether a vcov method exists for object x
    has_vcov_method <- paste0("vcov.", class(x)) %in% all_methods
    if (any(has_vcov_method)) {
      H <- -solve(vcov(x))
    } else {
      H <- NULL
    }
  }
  #
  # Set mle and nobs ----------
  #
  mle <- coef(x)
  n_obs <- nobs(x)
  #
  # Set V, using meat() or meatCL() from the sandwich package ----------
  #
  if (use_sandwich) {
    if (is.null(cluster)) {
      V <- sandwich::meat(x, fitted_object = x, loglik_fn = loglik_fn,
                          ...) * n_obs
    } else {
      V <- sandwich::meatCL(x, cluster = cluster, fitted_object = x,
                            loglik_fn = loglik_fn, ...) * n_obs
    }
  } else {
    V <- NULL
  }
  # We don't pass cluster because it would only be used in the estimation of
  # V: we have already estimated V using sandwich::meat() or sandwich::meatCL()
  result <- chandwich::adjust_loglik(loglik = loglik_fn,
                                     fitted_object = x,
                                     p = length(mle),
                                     par_names = names(mle),
                                     name = paste(class(x), collapse = "_"),
                                     mle = mle, H = H, V = V)
  class(result) <- c("oolax", "chandwich")
  return(result)
}
