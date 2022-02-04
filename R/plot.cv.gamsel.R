#' Plotting Routine for Gamsel Cross-Validation Object
#'
#' Produces a cross-validation curve with standard errors for a fitted gamsel
#' objecty.
#'
#' A plot showing cross-validation error is produced.  Nothing is returned.
#'
#' @param x \code{cv.gamsel} object
#' @param sign.lambda Either plot against \code{log(lambda)} (default) against
#' \code{-lambda} if \code{sign.lambda=-1}.
#' @param \dots Optional graphical parameters to plot.
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
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
#' gamsel.cv=cv.gamsel(X,y,bases=bases)
#' par(mfrow=c(1,1))
#' plot(gamsel.cv)
#' @method plot cv.gamsel
#' @export
plot.cv.gamsel=function(x,sign.lambda=1,...){
  cvobj=x
  xlab="log(Lambda)"
  if(sign.lambda<0)xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(cvobj$lambda),y=cvobj$cvm,ylim=range(cvobj$cvup,cvobj$cvlo),xlab=xlab,ylab=cvobj$name,type="n")
  new.args=list(...)
  if(length(new.args))plot.args[names(new.args)]=new.args
do.call("plot",plot.args)
error.bars(sign.lambda*log(cvobj$lambda),cvobj$cvup,cvobj$cvlo,width=0.01,col="darkgrey")
  points(sign.lambda*log(cvobj$lambda),cvobj$cvm,pch=20,col="red")
axis(side=3,at=sign.lambda*log(cvobj$lambda),labels=paste(cvobj$nz),tick=FALSE,line=0)
  lamin=c(cvobj$lambda.min,cvobj$lambda.1se)
  imin=c(cvobj$index.min,cvobj$index.1se)

  axis(side=3,at=sign.lambda*log(lamin),labels=paste(imin),tick=FALSE,line=-2,cex.axis=.7)
abline(v=sign.lambda*log(cvobj$lambda.min),lty=3)
abline(v=sign.lambda*log(cvobj$lambda.1se),lty=3)
  invisible()
}
