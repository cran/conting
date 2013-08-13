\name{plot.pval}
\alias{plot.pval}
\title{
Plot \code{pval} Objects
}
\description{
This function plots objects of class \code{"pval"}.
}
\usage{

\method{plot}{pval}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"pval"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}
\value{
This function will produce a plot of \code{T_obs} against \code{T_pred} (see \code{\link{bayespval}}) 
along with a line through the origin with slope one. The proportion of points above the line (grey points) 
gives the Bayesian p-value.
}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{For an example see \code{\link{bayespval}}.}

\seealso{
\code{\link{bayespval}}
}