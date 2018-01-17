\name{print.acceptrate}
\alias{print.acceptrate}

\title{
Prints \code{acceptrate} Objects
}
\description{
This function prints objects of class \code{"acceptrate"}.
}
\usage{

\method{print}{acceptrate}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class \code{"acceptrate"}.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
This function will simply print out the acceptance rates for the reversible jump and Metropolis-Hastings 
algorithms.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}

\note{For an example see \code{\link{accept_rate}}.}

\seealso{
\code{\link{accept_rate}}
}

