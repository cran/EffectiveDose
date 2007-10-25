\name{locpoly-class}
\alias{locpoly-class}
\alias{lines,locpoly-method}
\alias{plot,locpoly,missing-method}

\title{Class "locpoly" }
\description{ S4-Class for the value of the \code{locpoly} function}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("locpoly", ...)}. 
This class is the value of the function \code{locpoly} transfered to a S4-Class.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}. x-values of the smoothed fit in increasing order.}
    \item{\code{y}:}{Object of class \code{"numeric"}.The corresponding fitted y-values of the local polynomial estimator.}
    \item{\code{call}:}{Object of class \code{"call"}.  The call of the function \code{locpoly} for this object.} 
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "locpoly")}: line method to plot the class.}
    \item{plot}{\code{signature(x = "locpoly", y = "missing")}: to plot this class.}
  }
}


\seealso{
\code{\link[KernSmooth]{locpoly}}
}

\keyword{classes}
