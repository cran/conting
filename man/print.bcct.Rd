\name{print.bcct}
\alias{print.bcct}
\alias{print.bict}

\title{
Print \code{bcct} and \code{bict} Objects
}
\description{
This function prints objects of class \code{"bcct"} and \code{"bict"}.
}
\usage{

\method{print}{bcct}(x, ...)

\method{print}{bict}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{"bcct"} or \code{"bict"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
These functions print out very simple details on the \code{bcct} or \code{bict} objects. 
They display the number of cells in the table, the maximal model considered, the number of 
log-linear parameters in the maximal model, the number of MCMC iterations, the computer time 
required for the MCMC (in hours, minutes and seconds) and the prior used. 

In the case of objects of class \code{"bict"}, it also prints out the number of missing and 
censored cells.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk}.
}

\note{For examples see \code{\link{bcct}} and \code{\link{bict}}.}

\seealso{
\code{\link{bcct}},
\code{\link{bict}}.
}
