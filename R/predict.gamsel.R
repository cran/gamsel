#' Gamsel Prediction Routine
#'
#' Make predictions from a \code{gamsel} object.
#'
#'
#' @param object Fitted \code{gamsel} object.
#' @param newdata \code{nobs_new} x \code{p} matrix of new data values at which
#' to predict.
#' @param index Index of model in the sequence for which plotting is desired.
#' Note, this is NOT a lambda value.
#' @param type Type of prediction desired. Type \code{link} gives the linear
#' predictors for \code{"binomial"}, and fitted values for \code{"gaussian"}.
#' Type \code{response} gives fitted probabilities for \code{"binomial"} and
#' fitted values for \code{"gaussian"}. Type \code{"terms"} returns a matrix of
#' fitted functions, with as many columns as there are variables. Type
#' \code{nonzero} returns a list of the indices of nonzero coefficients at the
#' given \code{lambda} index.
#' @param \dots Not used
#' @return Either a vector aor a matrix is returned, depending on \code{type}.
#' @author Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
#' Hastie \email{hastie@@stanford.edu}
#' @seealso \code{\link{gamsel}}, \code{\link{cv.gamsel}},
#' \code{\link{summary.gamsel}}, \code{\link{basis.gen}}
#' @references Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
#' Additive Model Selection}
#' @method predict gamsel
#' @keywords regression smooth nonparametric
#' @examples
#'
#' ##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
#' data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
#' attach(data)
#' bases=pseudo.bases(X,degree=10,df=6)
#' # Gaussian gam
#' gamsel.out=gamsel(X,y,bases=bases)
#' preds=predict(gamsel.out,X,index=20,type="terms")
#'
#' @export
predict.gamsel=function(object, newdata, index=NULL,type=c("link","response","terms","nonzero"),...){
  type=match.arg(type)
  lambda=object$lambda
  nlambda=length(lambda)
  if(is.null(index))index=seq(nlambda) else  index=intersect(index,seq(nlambda))
  if(type=="nonzero")return(getActive(object,index,type="nonzero"))
  if(missing(newdata))stop("newdata is required for prediction")
  parms=object$parms
  betas=object$betas
  dimb=dim(betas)
  degrees=object$degrees
  p=length(degrees)
  dimx=dim(newdata)
  if(dimx[2]!=p)stop(paste("number of columns of x different from",p))
  U=as.list(1:p)
  offset=c(1,1+cumsum(degrees)[-p])
  for(i in 1:p)U[[i]]=basis.gen(newdata[,i],degree=degrees[i],parms=parms[[i]])
  betas=betas[,index,drop=FALSE]
  betas[offset,]=betas[offset,]+object$alpha[,index]
  dimnames(betas)=list(NULL,paste("l",index,sep=""))
    pred=switch(type,
        terms={  fitlist=as.list(1:p)
                 which=rep(1:p,degrees)
                 for(i in 1:p){
                    beta=betas[which==i,,drop=FALSE]
                    fitlist[[i]]=U[[i]]%*%beta
                  }
                 dd=c(dim(fitlist[[1]]),length(fitlist))
                 dn=c(dimnames(fitlist[[1]]),list(paste("v",1:p,sep="")))
                 array(do.call("cbind",fitlist),dd,dn)
              },
           {U=do.call("cbind",U)

              U%*%betas+rep(1,nrow(U)) %o% object$intercept[index]
            }
           )
  if(type=="response"&&object$family=="binomial"){
    pred=exp(pred)
    pred=pred/(1+pred)
  }
  pred
}
