% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.gen.R
\name{basis.gen}
\alias{basis.gen}
\title{Generate basis}
\usage{
basis.gen(x, df = 6, thresh = 0.01, degree = 8, parms = NULL, ...)
}
\arguments{
\item{x}{A vector of values for \code{basis.gen}, or a matrix for
\code{pseudo.bases}}

\item{df}{The degrees of freedom of the smoothing spline.}

\item{thresh}{If the next eigenvector improves the approximation by less
than threshold, a truncated bases is returned. For \code{pseudo.bases} this
can be a single value or a vector of values, which are recycled sequentially
for each column of \code{x}}

\item{degree}{The nominal number of basis elements. The basis returned has
no more than \code{degree} columns. For \code{pseudo.bases} this can be a
single value or a vector of values, which are recycled sequentially for each
column of \code{x}}

\item{parms}{A parameter set. If included in the call, these are used to
define the basis. This is used for prediction.}

\item{\dots}{other arguments}
}
\value{
the basis
}
\description{
Generate basis
}
