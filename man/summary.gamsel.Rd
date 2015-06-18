\name{summary.gamsel}
\alias{summary.gamsel}
\title{
 Gamsel summary routine
}
\description{
This makes a two-panel plot of the gamsel object.
}
\usage{
\method{summary}{gamsel}(object, label, ...)
}
\arguments{
  \item{object}{
\code{gamsel} object
}
  \item{label}{
if \code{TRUE}, annotate the plot with variable labels. Default is \code{FALSE}
}
  \item{\dots}{
additional arguments to summary
}
}
\details{
A two panel plot is produced, that summarizes the linear components and
the nonlinear components, as a function of lambda. For the linear
components, it is the coefficient for each variable. For the nonlinear,
we see the norm of the nonlinear coefficients.
}
\value{
  Nothing is returned.
}
\seealso{
\code{gamsel}, and methods \code{plot}, \code{print} and \code{predict}
for \code{cv.gamsel} object.
}
\references{
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized Additive Model
  Selection}
}
\author{Alexandra Chouldechova and Trevor Hastie\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\examples{
data=gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
attach(data)
bases=pseudo.bases(X,degree=10,df=6)
# Gaussian gam
gamsel.out=gamsel(X,y,bases=bases)
par(mfrow=c(1,2),mar=c(5,4,3,1))
summary(gamsel.out)
}
\keyword{regression}
\keyword{smooth}
\keyword{nonparametric}
