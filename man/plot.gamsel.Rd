% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.gamsel.R
\name{plot.gamsel}
\alias{plot.gamsel}
\title{Plotting Routine \code{gamsel} Object}
\usage{
\method{plot}{gamsel}(x, newx, index, which = 1:p, rugplot = TRUE, ylims, ...)
}
\arguments{
\item{x}{Fitted \code{gamsel} object.}

\item{newx}{\code{nobs_new} x \code{p} matrix giving values of each
predictor at which to plot.}

\item{index}{Index of lambda value (i.e., model) for which plotting is
desired.}

\item{which}{Which values to plot.  Default is all variables, i.e.
\code{{1,2,...,nvars}}. Besides indices, which can take two special values:
\code{"nonzero"} will plot only the nonzero functions, and
\code{"nonlinear"} only the nonlinear functions.}

\item{rugplot}{If \code{TRUE}, a rugplot showing values of \code{x} is shown
at the bottom of each fitted function plot.}

\item{ylims}{\code{ylim} argument for plotting each curve, which overides
the default which is the range of all the functions.}

\item{\dots}{Optional graphical parameters to plot.}
}
\description{
Produces plots of the estimated functions for specified variables at a given
value of \code{lambda}.
}
\details{
A plot of the specified fitted functions is produced.  Nothing is returned.
}
\examples{

##set.seed(1211)
##data=gamsel:::gendata(n=500,p=12,k.lin=3,k.nonlin=3,deg=8,sigma=0.5)
data = readRDS(system.file("extdata/gamsel_example.RDS", package = "gamsel"))
attach(data)
bases=pseudo.bases(X,degree=10,df=6)
# Gaussian gam
gamsel.out=gamsel(X,y,bases=bases)
par(mfrow=c(3,4))
plot(gamsel.out,newx=X,index=20)

}
\references{
Chouldechova, A. and Hastie, T. (2015) \emph{Generalized
Additive Model Selection}
}
\seealso{
\code{gamsel}, and \code{print.gamsel}, \code{summary.gamsel}
}
\author{
Alexandra Chouldechova and Trevor Hastie\cr Maintainer: Trevor
Hastie \email{hastie@stanford.edu}
}
\keyword{nonparametric}
\keyword{regression}
\keyword{smooth}
