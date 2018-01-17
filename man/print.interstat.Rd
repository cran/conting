\name{print.interstat}
\alias{print.interstat}

\title{
Print \code{interstat} Objects
}
\description{
This function prints objects of class \code{"interstat"}.
}
\usage{
\method{print}{interstat}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class \code{"interstat"}.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}
\value{
This function will print out the posterior probability, posterior mean, 
posterior variance and the 100*\code{prob.level}\% highest posterior density intervals (HPDIs) of each 
log-linear parameter (subject to the argument \code{cutoff}).
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\note{For an example see \code{\link{inter_stats}}.}

\seealso{
\code{\link{inter_stats}}.
}