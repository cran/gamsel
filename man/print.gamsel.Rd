\name{print.gamsel}
\alias{print.gamsel}
\title{
 print a gamsel object
}
\description{
Print a summary of the gamsel path at each step along the path
}
\usage{
\method{print}{gamsel}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted gamsel object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
  The call that produced the object \code{x} is printed, followed by a five-column
  matrix
  with columns \code{NonZero}, \code{Lin}, \code{NonLin}, \code{\%Dev}
  and \code{Lambda}. 
  The first three columns say how many nonzero, linear and nonlinear
  terms there are. \code{\%Dev} is the percent deviance
    explained (relative to the null deviance).
    }
\value{
The matrix above is silently returned}
\references{
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized Additive Model
  Selection}
}
\author{Alexandra Chouldechova and Trevor Hastie\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\seealso{
\code{\link{predict.gamsel}}, \code{\link{cv.gamsel}},
\code{\link{plot.gamsel}}, \code{\link{summary.gamsel}},
\code{\link{basis.gen}}, \code{\link{gendata}},  
}
\keyword{regression}
\keyword{smooth}
\keyword{nonparametric}
