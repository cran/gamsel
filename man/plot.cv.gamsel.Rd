\name{plot.cv.gamsel}
\alias{plot.cv.gamsel}
\title{
 Plotting Routine for Gamsel Cross-Validation Object
}
\description{
Produces a cross-validation curve with standard errors for a fitted
gamsel objecty.
}
\usage{
\method{plot}{cv.gamsel}(x, sign.lambda = 1, ...)
}
\arguments{
  \item{x}{
\code{cv.gamsel} object
}
  \item{sign.lambda}{
  Either plot against \code{log(lambda)} (default) against \code{-lambda} if \code{sign.lambda=-1}.
}
  \item{\dots}{
Optional graphical parameters to plot.
}
}
\details{
A plot showing cross-validation error is produced.  Nothing is returned.
}

\references{
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized Additive Model Selection}
}
\author{Alexandra Chouldechova and Trevor Hastie\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\examples{
data=gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
attach(data)
bases=pseudo.bases(X,degree=10,df=6)
# Gaussian gam
gamsel.out=gamsel(X,y,bases=bases)
gamsel.cv=cv.gamsel(X,y,bases=bases)
par(mfrow=c(1,1))
plot(gamsel.cv)
}
\keyword{regression}
\keyword{smooth}
\keyword{nonparametric}
