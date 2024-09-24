#' Fit Regularization Path for Generalized Additive Models
#'
#' Using overlap grouped lasso penalties, gamsel selects whether a
#' term in a gam is nonzero, linear, or a non-linear spline (up to a
#' specified max df per variable). It fits the entire regularization
#' path on a grid of values for the overall penalty lambda, both for
#' gaussian and binomial families. Key functions are [gamsel],
#' [predict.gamsel], [plot.gamsel], [print.gamsel], [summary.gamsel],
#' [cv.gamsel], [plot.cv.gamsel]
#'
#' @name gamsel-package
#' @useDynLib gamsel
#' @import splines
#' @import foreach
#' @import mda
#' @author Alexandra Chouldechova, Trevor Hastie Maintainer: Trevor Hastie
#' \email{hastie@@stanford.edu}
#' @keywords generalized additive models regression package
"_PACKAGE"





