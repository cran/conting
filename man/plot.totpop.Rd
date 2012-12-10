\name{plot.totpop}
\alias{plot.totpop}

\title{
Plot \code{totpop} Objects
}
\description{
This function plots objects of class "\code{totpop}".
}
\usage{
\method{plot}{totpop}(x, ...)
}

\arguments{
  \item{x}{
An object of class "\code{totpop}".
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
This function will produce a histogram of the MCMC sample from the posterior distribution of the 
total population size.
}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}

\note{For an example see \code{\link{total_pop}}.}

\seealso{
\code{\link{total_pop}}
}