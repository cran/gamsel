#' print a gamsel object
#'
#' Print a summary of the gamsel path at each step along the path
#'
#' The call that produced the object `x` is printed, followed by a
#' five-column matrix with columns `NonZero`, `Lin`, `NonLin`, `%Dev`
#' and `Lambda`.  The first three columns say how many nonzero, linear
#' and nonlinear terms there are. `%Dev` is the percent deviance
#' explained (relative to the null deviance).
#'
#' @param x fitted gamsel object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @return The matrix above is silently returned
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso [predict.gamsel], [cv.gamsel],
#' [plot.gamsel], [summary.gamsel],
#' [basis.gen]
#' @references Chouldechova, A. and Hastie, T. (2015) _Generalized Additive Model Selection_
#' @method print gamsel
#' @keywords regression smooth nonparametric
#' @export
print.gamsel=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      out=cbind(summarynz(x)[,c(3,1,2)],signif(x$dev.ratio,digits),signif(x$lambda,digits))
      colnames(out)=c("NonZero","Lin","NonLin","%Dev","Lambda")
      print(out,...)
      invisible(out)
    }
