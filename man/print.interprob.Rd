\name{print.interprob}
\alias{print.interprob}

\title{
Print \code{interprob} Objects
}
\description{
This function prints objects of class \code{"interprob"}.
}
\usage{
\method{print}{interprob}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class \code{"interprob"}.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
This function will print out the posterior probability of each term (subject 
to the argument \code{cutoff}).
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\note{For an example see \code{\link{inter_probs}}.}

\seealso{
\code{\link{inter_probs}}.
}
