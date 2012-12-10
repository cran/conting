\name{print.modprobs}
\alias{print.modprobs}

\title{
Print \code{modprobs} Objects
}
\description{
This function prints objects of class "\code{modprobs}".
}
\usage{
\method{print}{modprobs}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{
An object of class "\code{modprobs}".
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
This function will print out the posterior model probability of the 
"best" models as defined by the arguments \code{best} or \code{scale}.
}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{For an example see \code{\link{mod_probs}}.}

\seealso{
\code{\link{mod_probs}}.
}