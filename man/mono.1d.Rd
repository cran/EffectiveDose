\name{mono.1d}
\alias{mono.1d}
\title{strictly monotone and smooth function}
\description{
  this function applies a kernel smoothing method to monotonize a given fit with one independent variable
}
\usage{
mono.1d(fit, bandwidth, xx, kernel = "epanech", mono1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a list containing x-values and their corresponding y-values. The length of \code{fit[[1]]} determines degree of acuteness of the monotonizing procedure}
  \item{bandwidth}{a single number which represents the kernel bandwidth smoothing parameter. Missing values are not accepted.}
  \item{xx}{an additional vector of x-values where the monotonizing procedure is to be evaluated. If missing \code{fit[[1]]} is used instead.}
  \item{kernel}{ "`epanech"' - the Epanechnikov kernel}
  \item{mono1}{ either "increasing" or "decreasing"}
}


\value{
  returns an object of class "monofit"}
\references{Dette, H., Neumeyer, N., and Pilz, K. (2004) \emph{A simple nonparametric estimator of a monotone regression function}.
}

\author{ Regine Scheder \email{Regine.Scheder@rub.de}}

\examples{
data(cars)
speed<-cars$speed
dist<-cars$dist
fit1<-ksmooth(speed, dist, "normal", bandwidth=2)
##computes the Nadaraya-Watson estimate
fit2<-mono.1d(list(x=fit1$x, y=fit1$y),bandwidth=0.7, mono1="increasing") 				##calculates the monotone estimates
         plot(speed, dist)
         lines(fit1, col=2)
         lines(fit2, col=3)
}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}