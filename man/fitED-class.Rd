\name{fitED-class}
\docType{class}
\alias{fitED-class}
\alias{lines,fitED-method}
\alias{plot,fitED,missing-method}
\alias{print,fitED-method}
\alias{show,fitED-method}
\title{Class "fitED" }
\description{A object of class \code{fitED} is the ouput of \code{\link{ED}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("fitED", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Object of class \code{"numeric"}. Alpha-values where the effective dose levels is computed.}
    \item{\code{ED}:}{Object of class \code{"numeric"}. The levels of the Effective dose corresponding to the alpha-values.}
    \item{\code{call}:}{Object of class \code{"call"}.  The call of the function \code{"ED"} for this object.}
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "fitED")}: ... }
    \item{plot}{\code{signature(x = "fitED", y = "missing")}: ... }
    \item{print}{\code{signature(x = "fitED")}: ... } 
  }
}

\seealso{
  \code{\link[=ED.Boot.CI-class]{ED.Boot.CI}}, \code{\link{ED}}
}

\keyword{classes}
