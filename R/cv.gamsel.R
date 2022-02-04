#' Cross-validation Routine for Gamsel
#'
#' A routine for performing K-fold cross-validation for gamsel.
#'
#' This function has the effect of running \code{gamsel} \code{nfolds}+1 times.
#' The initial run uses all the data and gets the \code{lambda} sequence.  The
#' remaining runs fit the data with each of the folds omitted in turn.  The
#' error is accumulated, and the average error and standard deviation over the
#' folds is computed.  Note that \code{cv.gamsel} does NOT search for values
#' for \code{gamma}. A specific value should be supplied, else \code{gamma=.4}
#' is assumed by default. If users would like to cross-validate \code{gamma} as
#' well, they should call \code{cv.gamsel} with a pre-computed vector
#' \code{foldid}, and then use this same fold vector in separate calls to
#' \code{cv.gamsel} with different values of \code{gamma}. Note also that the
#' results of \code{cv.gamsel} are random, since the folds are selected at
#' random. Users can reduce this randomness by running \code{cv.gamsel} many
#' times, and averaging the error curves.
#'
#' @param x \code{x} matrix as in \code{gamsel}
#' @param y response \code{y} as in \code{gamsel}
#' @param lambda Optional use-supplied lambda sequence.  If \code{NULL},
#' default behaviour is for \code{gamsel} routine to automatically select a
#' good lambda sequence.
#' @param family \code{family} as in \code{gamsel}
#' @param degrees \code{degrees} as in \code{gamsel}
#' @param dfs \code{dfs} as in \code{gamsel}
#' @param bases \code{bases} as in \code{gamsel}
#' @param type.measure Loss function for cross-validated error calculation.
#' Currently there are four options: \code{mse} (mean squared error),
#' \code{mae} (mean absolute error), \code{deviance} (deviance, same as
#' \code{mse} for \code{family="gaussian"}), \code{class} (misclassification
#' error, for use with \code{family="binomial"}).
#' @param nfolds Numer of folds (default is 10).  Maximum value is \code{nobs}.
#' Small values of \code{nfolds} are recommended for large data sets.
#' @param foldid Optional vector of length \code{nobs} with values between 1
#' and \code{nfolds} specifying what fold each observation is in.
#' @param keep If \code{keep=TRUE}, a \emph{prevalidated} array is returned
#' containing fitted values for each observation and each value of
#' \code{lambda}. This means these fits are computed with this observation and
#' the rest of its fold omitted. The \code{folid} vector is also returned.
#' Default is keep=FALSE
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each
#' fold. See the example below for usage details.
#' @param \dots Other arguments that can be passed to \code{gamsel}.
#' @return an object of class \code{"cv.gamsel"} is returned, which is a list
#' with the ingredients of the cross-validation fit.  \item{lambda}{the values
#' of \code{lambda} used in the fits.} \item{cvm}{The mean cross-validated
#' error - a vector of length \code{length(lambda)}.} \item{cvsd}{estimate of
#' standard error of \code{cvm}.} \item{cvup}{upper curve = \code{cvm+cvsd}.}
#' \item{cvlo}{lower curve = \code{cvm-cvsd}.} \item{nzero}{number of non-zero
#' coefficients at each \code{lambda}.} \item{name}{a text string indicating
#' type of measure (for plotting purposes).} \item{gamsel.fit}{a fitted gamsel
#' object for the full data.} \item{lambda.min}{value of \code{lambda} that
#' gives minimum \code{cvm}.} \item{lambda.1se}{largest value of \code{lambda}
#' such that error is within 1 standard error of the minimum.}
#' \item{fit.preval}{if \code{keep=TRUE}, this is the array of prevalidated
#' fits. Some entries can be \code{NA}, if that and subsequent values of
#' \code{lambda} are not reached for that fold} \item{foldid}{if
#' \code{keep=TRUE}, the fold assignments used} \item{index.min}{the sequence
#' number of the minimum lambda.} \item{index.1se}{the sequence number of the
#' 1se lambda value.}
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso \code{gamsel}, \code{plot} function for \code{cv.gamsel} object.
#' @references Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
#' Additive Model Selection}
#' @keywords regression smooth nonparametric
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
#'
#' @export cv.gamsel
cv.gamsel <-
  function (x, y, lambda = NULL, family = c("gaussian", "binomial"),
            degrees = rep(10, p), dfs = rep(5, p), bases = pseudo.bases(x,
                                                     degrees, dfs,parallel=parallel, ...), type.measure = c("mse", "mae", "deviance",
                                                                           "class"), nfolds = 10, foldid, keep = FALSE, parallel = FALSE,
            ...)
{
  this.call = match.call()
  family = match.arg(family)
  p=ncol(x)
  y = drop(y)
  if (missing(type.measure))
    type.measure = switch(family, gaussian = "mse", binomial = "deviance")
  else type.measure = match.arg(type.measure)
  if (any(match(type.measure, c("deviance", "class"), FALSE)) &
      family == "gaussian") {
    warning(paste(type.measure, "not available for gaussian family; will use mse instead"))
    type.measure = "mse"
  }
  typenames = c(mse = "Mean-Squared Error", mae = "Mean Absolute Error",
    deviance = "Binomial Deviance", class = "Misclassification Error")
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.gamsel")
  nobs = nrow(x)
  y = drop(y)
  gamsel.object = gamsel(x, y, lambda = lambda, bases = bases, family=family,
    ...)
  lambda = gamsel.object$lambda
  if (missing(foldid))
    foldid = sample(rep(seq(nfolds), length = nobs))
  else nfolds = max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
#  if (parallel && require(foreach)) {
    outlist = foreach(i = seq(nfolds), .packages = c("gamsel")) %dopar%
    {
      which = foldid == i
      bases.sub = lapply(bases, basis.subset, subset = !which)
      gamsel(x[!which, , drop = FALSE], y[!which],
             lambda = lambda, bases = bases.sub, family=family,...)
    }
  }
  else {
    for (i in seq(nfolds)) {
      which = foldid == i
      bases.sub = lapply(bases, basis.subset, subset = !which)
      outlist[[i]] = gamsel(x[!which, , drop = FALSE],
               y[!which], lambda = lambda, bases = bases.sub, family=family,
               ...)
    }
  }
  predmat = matrix(NA, length(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, x[which, , drop = FALSE], type = "response")
    nlami = length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  prob_min = 1e-04
  prob_max = 1 - prob_min
  N = nobs - apply(is.na(predmat), 2, sum)
  cvraw = switch(type.measure, mse = (y - predmat)^2, mae = abs(y -
                                                        predmat), deviance = {
                                                          predmat = pmin(pmax(predmat, prob_min), prob_max)
                                                          lp = y * log(predmat) + (1 - y) * log(1-predmat)
                                                          -2 * lp
                                                        }, class = y * (predmat < 0.5) + (1 - y) * (predmat >= 0.5))
  nz = sapply(predict(gamsel.object, type = "nonzero"), length)
  cvm = apply(cvraw, 2, mean, na.rm = TRUE,trim=0.01)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE,trim=0.01)/(N -
                                                           1))
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
    cvsd, cvlo = cvm - cvsd, nzero = nz, name = typenames[type.measure],
    gamsel.fit = gamsel.object)
  if (keep)
    out = c(out, list(fit.preval = predmat, foldid = foldid))
  lamin = getmin(lambda, cvm, cvsd)
  imin = match(lamin, lambda)
  names(imin) = c("index.min", "index.1se")
  obj = c(out, as.list(c(lamin, imin)))
  class(obj) = "cv.gamsel"
  obj
}
