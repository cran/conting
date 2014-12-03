\name{print.totpop}
\alias{print.totpop}

\title{
Print \code{totpop} Objects
}
\description{
This function prints objects of class \code{"totpop"}.
}
\usage{
\method{print}{totpop}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class \code{"totpop"}.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
The function will print out the posterior mean and the 100*\code{prob.level}\% highest posterior 
density interval for the total population size.
}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk}.
}
\note{For an example see \code{\link{total_pop}}.}

\seealso{
\code{\link{total_pop}}.
}