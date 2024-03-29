% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.gamsel.R
\name{predict.gamsel}
\alias{predict.gamsel}
\title{Gamsel Prediction Routine}
\usage{
\method{predict}{gamsel}(
  object,
  newdata,
  index = NULL,
  type = c("link", "response", "terms", "nonzero"),
  ...
)
}
\arguments{
\item{object}{Fitted \code{gamsel} object.}

\item{newdata}{\code{nobs_new} x \code{p} matrix of new data values at which
to predict.}

\item{index}{Index of model in the sequence for which plotting is desired.
Note, this is NOT a lambda value.}

\item{type}{Type of prediction desired. Type \code{link} gives the linear
predictors for \code{"binomial"}, and fitted values for \code{"gaussian"}.
Type \code{response} gives fitted probabilities for \code{"binomial"} and
fitted values for \code{"gaussian"}. Type \code{"terms"} returns a matrix of
fitted functions, with as many columns as there are variables. Type
\code{nonzero} returns a list of the indices of nonzero coefficients at the
given \code{lambda} index.}

\item{\dots}{Not used}
}
\value{
Either a vector aor a matrix is returned, depending on \code{type}.
}
\description{
Make predictions from a \code{gamsel} object.
}
\examples{

##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
attach(data)
bases=pseudo.bases(X,degree=10,df=6)
# Gaussian gam
gamsel.out=gamsel(X,y,bases=bases)
preds=predict(gamsel.out,X,index=20,type="terms")

}
\references{
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
Additive Model Selection}
}
\seealso{
\code{\link{gamsel}}, \code{\link{cv.gamsel}},
\code{\link{summary.gamsel}}, \code{\link{basis.gen}}
}
\author{
Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
Hastie \email{hastie@stanford.edu}
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
