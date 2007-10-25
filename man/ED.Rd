\name{ED}
\docType{methods}
\alias{ED}
\alias{ED-methods}
\alias{ED,list,numeric-method}
\alias{ED,locfit,numeric-method}
\alias{ED,locpoly,numeric-method}
\alias{ED,list,missing-method}
\alias{ED,locfit,missing-method}
\alias{ED,locpoly,missing-method}


\title{Estimate for the Effective Dose level}
\description{
 estimates the effective dose level in binary response models
}
\usage{
ED(fitprob, alpha, ...)
}
\arguments{
  \item{fitprob}{ an object either of class \code{"list"}, \code{"locfit"}, or \code{"locpoly"}.}
  \item{alpha}{ alpha-levels where the effective dose is to be evaluated. Missing values are not accepted.}
  \item{\dots}{further parameters see details.}
}

\details{

The argument \code{"bandwidth"} is a single smoothing parameter of class \code{"numeric"} to evaluate the effective dose level. If missing a kind of Silverman's rule of thumb is used, i.e. \code{sd(x)/n^1/5}, where \code{x} is the dose level variable. 

The argument  \code{"N"} specifies the gridsize to approximate the integral to obtain the effective dose level. However, do not use too large values for \code{N}. By default \code{N=101} is used.

The argument  \code{"mono"} specifies if the effective dose level is assumed to be monotone increasing or decreasing. The default value is \code{mono="increasing"}, whereas in toxicology applications \code{mono="decreasing"} has to be used in order to get reasonable results.

The argument \code{"type"} can be used if the data consists of continuous measurements and the user still wants to the quantiles. Default is \code{type="cont"}. To obtain .5 quantile for continuous data use \code{type="prob"}.
}

\value{
 An object of class \code{\link[=fitED-class]{fitED}}.
 }
\references{ Dette, H., Neumeyer, N., and Pilz, K. (2005). \emph{A note on nonparametric estimation of the effective dose in quantal bioassay}. Journal of the American Statistical Association 100, 503-510. }
\author{ Regine Scheder \email{Regine.Scheder@rub.de}}

\seealso{ \code{\link[=fitED-class]{fitED}}. }
\examples{
#Function to generate data for a binary response model. The response variable has a success probability pnorm(x, mean=0.5, sd=0.5).
ybin=function(x){
	n=length(x)
	y=numeric(n)
	p=pnorm(x, mean=0.5, sd=0.5)
	for(i in 1:n){
	y[i]=rbinom(1,1,prob=p[i])
	}
	return(y)
}

x=seq(0,1,length.out=50) #the variable x presents the different dose levels
y=ybin(x)                              #y is the binary response simulated by the above function


fit=locfit(y~lp(x,deg=1, h=0.1, nn=0))   #fits a local linear estimate for the success probability using the locfit function
fit2=locpoly(x,y, degree=1, bandwidth=0.1) #fits a local linear estimate for the success probability using the locpoly function

#the effective dose level is estimated through 3 different ways.  
res=ED(list(x,y), alpha=seq(0.2,0.8, length.out=40)) #a list is used and internally the function locfit
res2=ED(fit, alpha=seq(0.2,0.8, length.out=40)) #the function ED is applied to a personally adjusted locfit object
res3=ED(fit2, alpha=seq(0.2,0.8, length.out=40)) #the function ED is applied to a personally adjusted locpoly object

#To display the results 
plot(res)
lines(res2, lty=2)
lines(res3,lty=3)
lines(pnorm(x, mean=0.5, sd=0.5), x, col="red") #adds the true effective dose levels 



}

\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}