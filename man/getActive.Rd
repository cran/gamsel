\name{getActive}
\alias{getActive}
\title{
 Returns active variables
}
\description{
Extract active variables of different kinds from a gamsel object
}
\usage{
getActive(object, index, type = , EPS = 0)
}
\arguments{
  \item{object}{
gamsel object
}
  \item{index}{
index or vector of indices at which to obtain active
information. \code{NULL} returns all.
}
  \item{type}{
type of active variables to report. One of \code{c("nonzero", "linear", "nonlinear")}
}
  \item{EPS}{
threshold for what is nonzero; default is 0
}
}
\details{
Returns a vector of variables indices of variables having the desired properties.
}
\value{
vector of indices
}