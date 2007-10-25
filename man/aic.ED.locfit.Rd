\name{aic.ED.locfit}
\docType{methods}
\alias{aic.ED.locfit}
\alias{aic.ED.locfit-methods}
\alias{aic.ED.locfit,fitED-method}

\title{Compute Akaike's Information Criterion.}
\description{
The function \code{aic.ED.locfit} returns the Akaike's information criterion for the effective dose level, where the \code{locfit} function is used for computing the estimate for the probability function.

The definition of AIC used here is -2*log-likelihood + pen*(fitted d.f.), where the fitted d.f. is used from the local linear fit, where the influence values define the degrees of freedom. To construct the log-likelihood basically the residual sum of squares is used.
}
\section{Methods}{
\describe{

\item{fit = "fitED"}{ value of the function \code{ED}, which is constructed using the \code{locfit} function. }
}}

\details{
The argument \code{"pen"} specifies the penalty for the degrees of freedom term. Default value is 2.

This function is geared to the function \code{aic} from the locfit package. Basically the AIC value is computed in the same way. To compute the log-likelihood value for the \code{ED} estimate the residual sum of squares is used, i.e. log-likelihood=-1/2*RSS as in the \code{locfit} package. The degrees of freedom are inherited from the locfit object. 
}

\value{
	\code{aic.ED.locfit} gives a numeric vector with 5 entries: \code{df} (degrees of freedom), \code{lk LL} (log-likelihood value of the locfit object), \code{lk ED} (likelihood value of the ED estimate), \code{AIC LL} (AIC score for the locfit object), \code{AIC ED} (AIC score for the ED estimate).
}

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

aic(fit)
res=ED(list(x,y), alpha=seq(0.2,0.8, length.out=40)) #a list is used and internally the function locfit
res2=ED(fit, alpha=seq(0.2,0.8, length.out=40)) #the function ED is applied to a personally adjusted \code{locfit} object

aic.ED.locfit(res)
aic.ED.locfit(res2)
}
\keyword{methods}
