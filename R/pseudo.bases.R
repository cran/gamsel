#' Generate pseudo-spline bases
#'
#' Generate an approximation to the Demmler-Reinsch orthonormal bases for
#' smoothing splines, using orthogonal polynomials. \code{basis.gen} generates
#' a basis for a single \code{x}, and \code{pseudo.bases} generates a list of
#' bases for each column of the matrix \code{x}.
#'
#' \code{basis.gen} starts with a basis of orthogonal polynomials of total
#' \code{degree}. These are each smoothed using a smoothing spline, which
#' allows for a one-step approximation to the Demmler-Reinsch basis for a
#' smoothing spline of rank equal to the degree. See the reference for details.
#' The function also approximates the appropriate diagonal penalty matrix for
#' this basis, so that the a approximate smoothing spline (generalized ridge
#' regression) has the target df.
#'
#' @param x A vector of values for \code{basis.gen}, or a matrix for
#' \code{pseudo.bases}
#' @param degree The nominal number of basis elements. The basis returned has
#' no more than \code{degree} columns. For \code{pseudo.bases} this can be a
#' single value or a vector of values, which are recycled sequentially for each
#' column of \code{x}
#' @param df The degrees of freedom of the smoothing spline.
#' @param parallel if TRUE, parallelize
#' @param \dots other arguments for \code{basis.gen} can be passed through to \code{[basis.gen]}
#' @return An orthonormal basis is returned (a list for \code{pseudo.bases}).
#' This has an attribute \code{parms}, which has elements
#' \code{coefs}Coefficients needed to generate the orthogonal polynomials
#' \code{rotate}Transformation matrix for transforming the polynomial basis
#' \code{d}penalty values for the diagonal penalty \code{df}df used
#' \code{degree}number of columns
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @references T. Hastie \emph{Pseudosplines}. (1996) JRSSB 58(2), 379-396.\cr
#' Chouldechova, A. and Hastie, T. (2015) \emph{Generalized Additive Model
#' Selection}
#' @keywords regression smooth nonparametric
#' @examples
#'
#' ##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
#' data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
#' attach(data)
#' bases=pseudo.bases(X,degree=10,df=6)
#' \dontrun{
#'      require(doMC)
#'      registerDoMC(cores=4)
#'      bases=pseudo.bases(X,degree=10,df=6,parallel=TRUE)
#' }
#'
#' @export pseudo.bases
pseudo.bases=function(x,degree=8,df=6,parallel=FALSE,...){
  p=ncol(x)
  degree=rep(degree,length=p)
  df=rep(df,length=p)
  out=as.list(1:p)

#  if (parallel && require(foreach)) {
  if (parallel){
    out = foreach(i = seq(p), .packages = c("gamsel")) %dopar%
    {
      basis.gen(x[,i],degree=degree[i],df=df[i],...)
    }
  }
  else {
    for(i in 1:p){
      out[[i]]=basis.gen(x[,i],degree=degree[i],df=df[i],...)
    }
  }
  out
}

