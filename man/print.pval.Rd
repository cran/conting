\name{print.pval}
\alias{print.pval}

\title{
Print \code{pval} Objects
}
\description{
This function prints objects of class \code{"pval"}.
}
\usage{
\method{print}{pval}(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"pval"}.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
This function will print out summary statistics for the predictive (\code{T_pred}) and 
observed (\code{T_obs}) discrepancy statistics. Additionally it will output the associated 
Bayesian p-value.
}

\author{
Antony M. Overstall (\email{antony@mcs.st-and.ac.uk}).
}

\note{For an example see \code{\link{bayespval}}.}

\seealso{
\code{\link{bayespval}}
}
