#' Fit Regularization Path for Gaussian or Binomial Generalized Additive Model
#'
#' Using overlap grouped lasso penalties, gamsel selects whether a term in a
#' gam is nonzero, linear, or a non-linear spline (up to a specified max df per
#' variable). It fits the entire regularization path on a grid of values for
#' the overall penalty lambda, both for gaussian and binomial families.
#'
#' The sequence of models along the \code{lambda} path is fit by (block)
#' cordinate descent.  In the case of logistic regression the fitting routine
#' may terminate before all \code{num_lambda} values of \code{lambda} have been
#' used.  This occurs when the fraction of null deviance explained by the model
#' gets too close to 1, at which point the fit becomes numerically unstable.
#' Each of the smooth terms is computed using an approximation to the
#' Demmler-Reinsch smoothing spline basis for that variable, and the
#' accompanying diagonal pernalty matrix.
#'
#' @param x Input (predictor) matrix of dimension \code{nobs} x \code{nvars}.
#' Each observation is a row.
#' @param y Response variable.  Quantitative for \code{family="gaussian"} and
#' with values in \code{{0,1}} for \code{family="binomial"}
#' @param num_lambda Number of \code{lambda} values to use. (Length of
#' \code{lambda} sequence.)
#' @param lambda User-supplied \code{lambda} sequence.  For best performance,
#' leave as \code{NULL} and allow the routine to automatically select
#' \code{lambda}.  Otherwise, supply a (preferably gradually) decreasing
#' sequence.
#' @param family Response type. \code{"gaussian"} for linear model (default).
#' \code{"binomial"} for logistic model.
#' @param degrees An integer vector of length \code{nvars} specifying the
#' maximum number of spline basis functions to use for each variable.
#' @param gamma Penalty mixing parameter \eqn{0 \le\gamma\le 1}.  Values \eqn{
#' \gamma < 0.5} penalize linear fit less than non-linear fit.  The default is
#' \eqn{\gamma = 0.4}, which encourages a linear term over a nonlinear term.
#' @param dfs Numeric vector of length \code{nvars} specifying the maximum
#' (end-of-path) degrees of freedom for each variable.
#' @param bases A list of orthonormal bases for the non-linear terms for each
#' variable. The function \code{pseudo.bases} generates these, using the
#' parameters \code{dfs} and \code{degrees}. See the documentation for
#' \code{\link{pseudo.bases}}.
#' @param tol Convergence threshold for coordinate descent.  The coordinate
#' descent loop continues until the total change in objective after a pass over
#' all variables is less than \code{tol}.  Default is \code{1e-4}.
#' @param max_iter Maximum number of coordinate descent iterations over all the
#' variables for each \code{lambda} value.  Default is 2000.
#' @param traceit If \code{TRUE}, various information is printed during the
#' fitting process.
#' @param parallel passed on to the \code{pseudo.bases()} function. Uses
#' multiple process if available.
#' @param \dots additional arguments passed on to \code{pseudo.bases()}
#' @return An object with S3 class \code{gamsel}.  %% If it is a LIST, use
#' \item{intercept}{Intercept sequence of length \code{num_lambda}}
#' \item{alphas}{\code{nvars} x \code{num_lambda} matrix of linear coefficient
#' estimates} \item{betas}{\code{sum(degrees)} x \code{num_lambda} matrix of
#' non-linear coefficient estimates} \item{lambdas}{The sequence of lambda
#' values used} \item{degrees}{Number of basis functions used for each
#' variable} \item{parms}{A set of parameters that capture the bases used. This
#' allows for efficient generation of the bases elements for
#' \code{predict.gamsel}}, the \code{predict} method for this class.
#' \item{family}{\code{"gaussian"} or \code{"binomial"}} \item{nulldev}{Null
#' deviance (deviance of the intercept model)} \item{dev.ratio}{Vector of
#' length \code{num_lambda} giving fraction of (null) deviance explained by
#' each model along the \code{lambda} sequence} \item{call}{The call that
#' produced this object} %% ...
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso \code{\link{predict.gamsel}}, \code{\link{cv.gamsel}},
#' \code{\link{plot.gamsel}}, \code{\link{summary.gamsel}},
#' \code{\link{basis.gen}},
#' @references Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
#' Additive Model Selection}, \url{https://arxiv.org/abs/1506.03850}
#' @keywords regression smooth nonparametric
#' @importFrom graphics abline axis lines plot points rug segments text
#' @importFrom grDevices rainbow
#' @importFrom utils head tail packageDescription
#' @importFrom stats poly polym predict rbinom rnorm
#'
#' @examples
#'
#' ##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
#' data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
#' attach(data)
#' bases=pseudo.bases(X,degree=10,df=6)
#' # Gaussian gam
#' gamsel.out=gamsel(X,y,bases=bases)
#' par(mfrow=c(1,2),mar=c(5,4,3,1))
#' summary(gamsel.out)
#' gamsel.cv=cv.gamsel(X,y,bases=bases)
#' par(mfrow=c(1,1))
#' plot(gamsel.cv)
#' par(mfrow=c(3,4))
#' plot(gamsel.out,newx=X,index=20)
#' # Binomial model
#' gamsel.out=gamsel(X,yb,family="binomial")
#' par(mfrow=c(1,2),mar=c(5,4,3,1))
#' summary(gamsel.out)
#' par(mfrow=c(3,4))
#' plot(gamsel.out,newx=X,index=30)
#'
#' @export gamsel
gamsel <-
  function (x,y, num_lambda = 50,lambda=NULL, family = c("gaussian","binomial"), degrees = rep(10,p), gamma = 0.4, dfs=rep(5,p), bases=pseudo.bases(x,degrees,dfs,parallel=parallel,...), tol = 1e-04, max_iter = 2000,traceit=FALSE, parallel=FALSE, ...)
{
  this.call=match.call()
  family=match.arg(family)
  family_id <- as.integer(ifelse(family == "gaussian", 0, 1))
  n <- length(y)
  p <- ncol(x)
  ##Create U, X, D and psis  from the bases
  degrees=sapply(bases,dim)[2,]
  U=do.call("cbind",bases)
  X=do.call("cbind",lapply(bases,function(x)x[,1,drop=FALSE]))
  parms=lapply(bases,"attr","parms")
  getdpsi=function(parms){
    d=parms$d
    if(length(d)>1){
      psi=d[2]
      d=d/d[2]
      d[1]=1
    }
    else{
      d=1
      psi=0
    }
    list(d=d,psi=psi)
  }
  dpsi=lapply(parms,getdpsi)
  D_mat=unlist(lapply(dpsi,"[[","d"))
  psi=sapply(dpsi,"[[","psi")
  if(is.null(lambda)) lambda=rep(-0.5,num_lambda) else {
    lambda=as.numeric(lambda)
    lambda=rev(unique(sort(lambda)))
  }
  num_lambda <- as.integer(length(lambda))
  degrees <- as.integer(degrees)
  max_iter <- as.integer(max_iter)
  out <- .Call("gamselFit", y, X, U, tol, degrees, D_mat, gamma,
               psi, family_id, max_iter, lambda,num_lambda,as.integer(traceit))[c("intercept","alphas","betas","lambdas")]
  out$degrees=degrees
  out$parms=parms
  out$family=family
  out=c(out,fracdev(U,y,out$alphas,out$betas,out$intercept,degrees,family))
  out$call=this.call
  class(out)="gamsel"
  out
}
