\name{Boot.CI-methods}
\docType{methods}
\alias{Boot.CI}
\alias{Boot.CI-methods}
\alias{Boot.CI,list,numeric-method}
\alias{Boot.CI,locfit,numeric-method}
\alias{Boot.CI,locpoly,numeric-method}

\alias{Boot.CI,list,missing-method}
\alias{Boot.CI,locfit,missing-method}
\alias{Boot.CI,locpoly,missing-method}

\title{Bootstrap confidence intervals for ED levels}
\description{
 gives the Bootstrap confidence interval for the effective dose level.
}
\usage{
Boot.CI(fitprob, alpha, ...)
}
\arguments{
  \item{fitprob}{ an object either of class \code{"list"}, \code{"locfit"}, or \code{"locpoly"}.}
  \item{alpha}{alpha-levels where the bootstrap confidence intervals for  effecitve dose levels is to be estimated. Missing values are not accepted.}
  \item{\dots}{ further parameters see details. }
}
\details{
The confidence level of the bootstrap interval can be fixed through the argument \code{"level"}. The default is \code{level=0.95}.

The argument \code{"bandwidth"} is a single smoothing parameter of class \code{"numeric"} to evaluate the effective dose level. If missing a kind of Silverman's rule of thumb is used.

With the argument \code{"R"}, the number of bootstrap replications can be determined. By default, the procedure uses \code{R=100}.

 The argument  \code{"N"} specifies the gridsize to approximate the integral to obtain the effective dose level. However, do not use too large values for \code{N}. By default \code{N=101} is used.
 
 The argument  \code{"mono"}  specifies if the effective dose level is assumed to be monotone increasing or decreasing. The default value is \code{mono="increasing"}, whereas in toxicology applications \code{mono="decreasing"} has to be used in order to get reasonable results.

The argument \code{"type"} can be used if the data consists of continuous measurements and the user still wants to the quantiles. Default is \code{type="cont"}. To obtain .5 quantile for continuous data use \code{type="prob"}.

 Bootstrap confidence intervals might give inaccurate lower and upper bounds of the confidence interval, if the underlying dose response curve \code{"p"} is rather flat. The function \code{Boot.CI} gives a warning if the estimated derivative of the dose response curve is smaller than 0.1 in absolute values. 
}
\value{
  An object of class  \code{\link[=ED.Boot.CI-class]{ED.Boot.CI}}.
  }
\author{Regine Scheder \email{Regine.Scheder@rub.de}}


\seealso{ \code{\link{ED}}, \code{\link[=ED.Boot.CI-class]{ED.Boot.CI}}. }
\examples{
#Same model as discussed for the function ED
ybin=function(x){
	n=length(x)
	y=numeric(n)
	p=pnorm(x, mean=0.5, sd=0.5)
	for(i in 1:n){
	y[i]=rbinom(1,1,prob=p[i])
	}
	return(y)
}

x=seq(0,1,length.out=50)
y=ybin(x)

fit=locfit(y~lp(x,deg=1, h=0.1, nn=0))
fit2=locpoly(x,y, degree=1, bandwidth=0.1)

#For each usage of the function Boot.CI a Bootstrap confidence interval is computed

Boot.CI(list(x,y), alpha=0.5)
Boot.CI(fit,alpha=0.5)
Boot.CI(fit2,alpha=0.5)

#In this example, it is not useful to compute the Bootstrap confidence intervals close to the boundaries, since the function p is very flat
ybin=function(x){
        n=length(x)
        y=numeric(n)
        p=pnorm(x, mean=0.5, sd=0.05)
        for(i in 1:n){
        y[i]=rbinom(1,1,prob=p[i])
        }
        return(y)
}
#The following plot shows how flat the function p is at the bounds.
plot(x,pnorm(x, mean=0.5, sd=0.05), type="l")

x=seq(0,1,length.out=50)
y=ybin(x)

fit=locfit(y~lp(x,deg=1, h=0.1, nn=0))
fit2=locpoly(x,y, degree=1, bandwidth=0.1)

#For each usage of the function Boot.CI a Bootstrap confidence interval is computed and the upper and lower limits are returned followed by a warning.

Boot.CI(list(x,y), alpha=0.2)
Boot.CI(fit,alpha=0.2)
Boot.CI(fit2,alpha=0.2)

}

\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}