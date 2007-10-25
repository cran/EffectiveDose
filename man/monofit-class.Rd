\name{monofit-class}
\alias{monofit-class}
\alias{lines,monofit-method}
\alias{plot,monofit,missing-method}

\title{Class "monofit" }
\description{A object of class \code{monofit} is the ouput of \code{\link{mono.1d}}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("monofit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{x-variable is of class \code{"numeric"}}
    \item{\code{y}:}{y-variable is of class \code{"numeric"}}
    \item{\code{z}:}{z-variable is of class \code{"matrix"}}
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "monofit")}: ... }
    \item{plot}{\code{signature(x = "monofit", y = "missing")}: ... }
  }
}


\seealso{
  \code{\link{mono.1d}}
}

\keyword{classes}
