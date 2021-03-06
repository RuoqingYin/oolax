#' oolax: Objected-orientated Loglikelihood Adjustment for Extreme Value Models
#'
#'Provides functions for loglikelihood adjustment for fitting extreme value
#'models. The fitted object is from the package 'evd' and 'ismev' for extreme value
#'analysis. The 'chandwich' package
#'<https://cran.r-project.org/web/packages/chandwich/index.html> is used to construct
#'the adjustment given the clustered data. The 'sandwich' package
#'<https://cran.r-project.org/web/packages/sandwich/index.html> is used to adjust the
#'robust covariance matrix estimators.
#'
#'
#'
#' @details
#' The main function in the oolax package is \code{\link{alogLik}}, which provides
#' loglikelihood adjustment for fitting Generalised Extreme Value (GEV) model and
#' Threshold Modelling using generalised Pareto distribution (GPD). The adjustment is
#' based on package \code{\link[evd]{evd}} and \code{\link[ismev]{ismev}}. The user can
#' make adjustment towards fitted obejcts such as \link[evd]{fgev} and \link[evd]{fgev}
#' from the package \code{\link[evd]{evd}}, and \link[ismev]{gev.fit} and
#' \link[ismev]{gpd.fit} from the package \code{\link[ismev]{ismev}}. The metholdology is
#' develped based on \code{\link{sandwich}} and \code{\link{chandwich}}.
#'
#' See \code{vignette("oolax-vignette", package = "oolax")} for an overview and detailed
#' examples about the package.
#'
#'
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @docType package
#' @name oolax
#' @import methods
#' @import sandwich
NULL
