#' Gamsel summary routine
#'
#' This makes a two-panel plot of the gamsel object.
#'
#' A two panel plot is produced, that summarizes the linear components and the
#' nonlinear components, as a function of lambda. For the linear components, it
#' is the coefficient for each variable. For the nonlinear, we see the norm of
#' the nonlinear coefficients.
#'
#' @param object \code{gamsel} object
#' @param label if \code{TRUE}, annotate the plot with variable labels. Default
#' is \code{FALSE}
#' @param \dots additional arguments to summary
#' @return Nothing is returned.
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso \code{gamsel}, and methods \code{plot}, \code{print} and
#' \code{predict} for \code{cv.gamsel} object.
#' @references Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
#' Additive Model Selection}
#' @method summary gamsel
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
#'
#' @export
summary.gamsel <-
  function (object, label=FALSE,...)
{
  alphas=object$alphas
  betas=object$betas
  p=dim(alphas)[1]
  degrees=object$degrees
  maxvars = p
  col.end = 0.7
  colours = c(rainbow(n = maxvars, end = col.end, v = 0.9),
    rainbow(n = (p - maxvars), start = col.end + 0.1, s = 0.25))
  lambda=object$lambda

  lambda.max = max(lambda)
  lambda.min = min(lambda[lambda>0.001])
  nzlines=function(lambda,alpha,...){
    if(any(abs(alpha)>0)){
      num_lambda <- length(lambda)
      start=max(1,min(seq(num_lambda)[abs(alpha)>0])-1)
      whichnz=seq(from=start,to=num_lambda)
      if (length(whichnz)>1){
        lines(lambda[whichnz], alpha[whichnz],...)
      }
    }
    invisible()
  }
  plot(0, type = "n", xlab = expression(lambda), ylab = expression(alpha),
       xlim = c(lambda.max * 1.05, lambda.min * (0.9-.1*label)), ylim = range(alphas), main = "Linear Components",  log = "x")
  abline(h=0,lty=3)
  for (j in 1:maxvars) {
    nzlines(lambda,alphas[j,], col = colours[j], lwd = 2, type = "l", pch = "", cex = 0)
  }
  if(label)text(rep(lambda.min*.85,maxvars),alphas[,length(lambda)],labels=seq(maxvars),col=colours,cex=0.6)
  nbetas=dim(betas)[1]
  counter=rep(1:maxvars,degrees)
  counter=diag(maxvars)[counter,]
  nbeta=sqrt(t(counter)%*%(betas^2))
  betamax=max(nbeta)
  plot(0, type = "n", xlab = expression(lambda), ylab = expression(paste("||",beta,"||")),
       xlim = c(lambda.max * 1.05, lambda.min * (0.9-.1*label)), ylim = c(0,1.05 * betamax), main = "Non-linear Components",  log = "x")
  abline(h=0,lty=3)
  for (j in 1:maxvars) {
    nzlines(lambda,nbeta[j,], col = colours[j], lwd = 2, type = "l", pch = "")

  }
  if(label)text(rep(lambda.min*.85,maxvars),nbeta[,length(lambda)],labels=seq(maxvars),col=colours,cex=0.6)

}
