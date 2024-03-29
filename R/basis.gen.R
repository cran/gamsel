#' Generate basis
#' @inheritParams pseudo.bases
#' @param thresh If the next eigenvector improves the approximation by less
#' than threshold, a truncated bases is returned. For \code{pseudo.bases} this
#' can be a single value or a vector of values, which are recycled sequentially
#' for each column of \code{x}
#' @param parms A parameter set. If included in the call, these are used to
#' define the basis. This is used for prediction.
#' @param \dots other arguments
#' @return the basis
#' @export basis.gen
basis.gen <- function(x, df = 6, thresh=0.01, degree = 8, parms=NULL,...)
{
##uses the function mspline in the mda package, which fits a vector spline
    if(!is.null(parms)){#prediction
    poly(x,degree=parms$degree,coefs=parms$coefs)%*%parms$rotate
  } else
  {
    norm <- 1
    sdf <- df*1.5
    ratio <- NULL
    if(degree<1){
      warning("degree <1 converted to 1")
      degree=1
    }
    px <- poly(x, degree=degree)
    if(degree > 1){
      spx <- px
#      spx[,2:degree]=mspline(x,px[,2:degree,drop=FALSE],df=sdf)
      spx[,2:degree]=mspline(x,px[,2:degree,drop=FALSE],df=sdf)$s
      psp <- matrix(0, degree, degree)
      psp[1, 1] <- 1
      for(k in 2:degree) {
        for(i in 2:k)
          psp[i, k] <- min(1, sum(spx[, k] * px[, i]))
        normk <- 2 * sum(psp[1:k - 1, k]^2) + psp[k, k]^2
        norm <- c(norm, normk)
        wdof <- sum(norm)
        ratio <- norm[k]/wdof
        if(ratio < thresh)
          break
      }
      ndegree=length(norm)
      psp <- psp[seq(ndegree), seq(ndegree)]
      psp <- psp + t(psp) - diag(diag(psp))
      e.psp <- eigen(psp)
      d=(1/e.psp$values) -1
      lambdadf=df.inv(d,df,1)
      sbasis = px[, seq(ndegree)] %*% e.psp$vectors
      bcoefs=attr(px,"coefs")
      if(ndegree<degree){
        bcoefs$alpha=bcoefs$alpha[seq(ndegree)]
        bcoefs$norm2=bcoefs$norm2[seq(ndegree+2)]
      }
      parms=list(coefs=bcoefs,rotate=e.psp$vectors,d=d*lambdadf$lambda, df=lambdadf$df,degree=ndegree)
    }else
    { #degree =1
      sbasis=px
      parms=list(coefs=attr(px,"coefs"),rotate=matrix(1),d=1, df=1,degree=1)
    }
    attr(sbasis,"parms")=parms
    sbasis
  }
}



