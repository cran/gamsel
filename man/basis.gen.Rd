\name{basis.gen}
\alias{basis.gen}
\alias{pseudo.bases}
\title{Generate pseudo-spline bases}
\description{
Generate an approximation to the Demmler-Reinsch orthonormal bases for
smoothing splines, using orthogonal polynomials. \code{basis.gen}
generates a basis for a single \code{x}, and \code{pseudo.bases}
generates a list of bases for each column of the matrix \code{x}.
}
\usage{
basis.gen(x, df = 6, thresh = 0.01, degree = 8, parms = NULL,...)
pseudo.bases(x, degree = 8, df = 6, parallel=FALSE, ...)
}
\arguments{
  \item{x}{
A vector of values for \code{basis.gen}, or a matrix for \code{pseudo.bases}
}
  \item{df}{
The degrees of freedom of the smoothing spline.
}
  \item{thresh}{
If the next eigenvector improves the approximation by less than
threshold, a truncated bases is returned. For \code{pseudo.bases} this can be a single
value or a vector of values, which are recycled sequentially for each
column of \code{x}
}
  \item{degree}{
The nominal number of basis elements. The basis returned has no more
than \code{degree} columns. For \code{pseudo.bases} this can be a single
value or a vector of values, which are recycled sequentially for each
column of \code{x}
}
  \item{parms}{
A parameter set. If included in the call, these are used to define the
basis. This is used for prediction.
}
\item{parallel}{For \code{pseudo.bases}, allows for parallel bases
  computation in multiple cores.
}
\item{\dots}{other arguments for \code{basis.gen} can be passed through \code{pseudo.bases}}
}
\details{
\code{basis.gen} starts with a basis of orthogonal polynomials of total
\code{degree}. These are each smoothed using a smoothing spline, which
allows for a one-step approximation to the  Demmler-Reinsch basis for a
smoothing spline of rank equal to the degree. See the reference for
details. The function also
approximates the appropriate diagonal penalty matrix for this basis, so
that the a approximate smoothing spline (generalized ridge regression)
has the target df.
}
\value{
  An orthonormal basis is returned (a list for \code{pseudo.bases}).
  This has an attribute \code{parms}, which has elements
  \code{coefs}{Coefficients needed to generate the orthogonal
    polynomials}
  \code{rotate}{Transformation matrix for transforming the polynomial
    basis}
  \code{d}{penalty values for the diagonal penalty}
  \code{df}{df used}
  \code{degree}{number of columns}
}
\references{T. Hastie \emph{Pseudosplines}. (1996) JRSSB 58(2), 379-396.\cr
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized Additive Model Selection}
}
\author{Alexandra Chouldechova and Trevor Hastie\cr
Maintainer: Trevor Hastie \email{hastie@stanford.edu}}
\examples{
data=gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
attach(data)
bases=pseudo.bases(X,degree=10,df=6)
\dontrun{
     require(doMC)
     registerDoMC(cores=4)
     bases=pseudo.bases(X,degree=10,df=6,parallel=TRUE)
}
}
\keyword{regression}
\keyword{smooth}
\keyword{nonparametric}
