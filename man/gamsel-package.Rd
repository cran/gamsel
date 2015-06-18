\name{gamsel-package}
\alias{gamsel-package}
%\alias{gamsel}
\docType{package}
\title{gamsel}
\description{
Using overlap grouped lasso penalties, gamsel selects whether a term in a gam is nonzero, linear, or a non-linear spline (up to a specified max df per variable). It fits the entire regularization path on a grid of values for the overall penalty lambda, both for gaussian and binomial families. 
}
\details{
\tabular{ll}{
Package: \tab gamsel\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-05\cr
License: \tab What license is it under?\cr
}
Accepts \code{x, y} data and complexity/tuning parameters.  Key functions:
\code{gamsel}\cr
\code{predict.gamsel}\cr
\code{plot.gamsel}\cr
\code{print.gamsel}\cr
\code{summary.gamsel}\cr
\code{cv.gamsel}\cr
\code{plot.cv.gamsel}
}
\author{
Alexandra Chouldechova, Trevor Hastie
Maintainer:  Trevor Hastie \email{hastie@stanford.edu}
}
\keyword{generalized additive models}
\keyword{regression}
\keyword{package}
