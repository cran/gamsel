#' Plotting Routine \code{gamsel} Object
#'
#' Produces plots of the estimated functions for specified variables at a given
#' value of \code{lambda}.
#'
#' A plot of the specified fitted functions is produced.  Nothing is returned.
#'
#' @param x Fitted \code{gamsel} object.
#' @param newx \code{nobs_new} x \code{p} matrix giving values of each
#' predictor at which to plot.
#' @param index Index of lambda value (i.e., model) for which plotting is
#' desired.
#' @param which Which values to plot.  Default is all variables, i.e.
#' \code{{1,2,...,nvars}}. Besides indices, which can take two special values:
#' \code{"nonzero"} will plot only the nonzero functions, and
#' \code{"nonlinear"} only the nonlinear functions.
#' @param rugplot If \code{TRUE}, a rugplot showing values of \code{x} is shown
#' at the bottom of each fitted function plot.
#' @param ylims \code{ylim} argument for plotting each curve, which overides
#' the default which is the range of all the functions.
#' @param \dots Optional graphical parameters to plot.
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso \code{gamsel}, and \code{print.gamsel}, \code{summary.gamsel}
#' @references Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
#' Additive Model Selection}
#' @method plot gamsel
#' @keywords regression smooth nonparametric
#' @examples
#'
#' ##set.seed(1211)
#' ##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
#' data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
#' attach(data)
#' bases=pseudo.bases(X,degree=10,df=6)
#' # Gaussian gam
#' gamsel.out=gamsel(X,y,bases=bases)
#' par(mfrow=c(3,4))
#' plot(gamsel.out,newx=X,index=20)
#'
#' @export
plot.gamsel <-
  function(x,newx, index,which=1:p,rugplot=TRUE,ylims,...){
    gamsel.out=x
    x=newx
    if(missing(index))index=length(gamsel.out$lambdas)
    pJ=dim(gamsel.out$alpha)
    p=pJ[1]
    maxJ=pJ[2]
    nonlin <- getActive(gamsel.out, index,type="nonlinear")[[1]]
    active_linear <- getActive(gamsel.out, index,type="linear")[[1]]
    lin=setdiff(active_linear,nonlin)
    zero=setdiff(1:p,union(nonlin,lin))
    if(which[1]=="nonzero")which=seq(p)[-zero]
    if(which[1]=="nonlinear")which=nonlin
    if(is.null(which)){warning("Nothing to plot with this choice of argument which");return()}
    colvals=rep("blue",p)
    colvals[lin]="green"
    colvals[nonlin]="red"
    termmat=drop(predict.gamsel(gamsel.out,x,index=index,type="terms"))
    lambda=gamsel.out$lambda[index]
    if(missing(ylims))ylims=range(termmat[,which])
    for(ell in which){
      o=order(x[,ell])
      plot(x[o,ell], termmat[o,ell], type='l', ylab=paste("f(v",ell,")",sep=""),xlab=paste("v",ell,sep=""),ylim=ylims,col=colvals[ell],lwd=2,...)
      if(rugplot)rug(x[,ell])
    }
  }
