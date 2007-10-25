\name{ED.Boot.CI-class}
\alias{ED.Boot.CI-class}
\alias{print,ED.Boot.CI-method}
\alias{show,ED.Boot.CI-method}
\title{Class "ED.Boot.CI" }
\description{A object of class \code{ED.Boot.CI} is the ouput of \code{\link{Boot.CI}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ED.Boot.CI", ...)}.
}
\section{Slots}{
  \describe{
   \item{\code{CI}:}{Object of class \code{"numeric"}. Lower and upper limits of the bootstrap interval  labeled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5 \% and 97.5 \% ).}
    \item{\code{R}:}{Object of class \code{"numeric"}. Number of Bootstrap replications to calculate the confidence interval.}
    \item{\code{call}:}{Object of class \code{"call"}. The call of the function \code{"ED"} for this object.}
   }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "ED.Boot.CI")}: ... }
  }
}

\seealso{
  \code{\link[=fitED-class]{fitED}}, \code{\link{Boot.CI}}
}

\keyword{classes}
