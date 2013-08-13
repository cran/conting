\name{print.submod}
\alias{print.submod}

\title{
Print \code{submod} Objects
}
\description{
This function prints objects of class \code{"submod"}.
}
\usage{
\method{print}{submod}(x, ..., digits = max(3, getOption("digits") - 3))
}

\arguments{
  \item{x}{
An object of class \code{"submod"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
  \item{digits}{
An optional argument controlling the rounding of output.
}
}

\value{
Firstly, conditional on the model of interest (defined by \code{formula} and \code{order}), 
this function will print out the posterior means, posterior variances and 100*\code{prob.level}\% 
highest posterior density intervals (HPDIs) for each of the log-linear parameters. Secondly, conditional 
on the model of interest, it will print out summaries of the discrepancy statistics and the 
corresponding Bayesian p-value. Finally, if the class of the object passed to \code{sub_model} is 
\code{"bict"}, then it will print out the posterior mean and 100*\code{prob.level}\% HPDI for the total
population size, conditional on the model of interest.

}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{For an example see \code{\link{sub_model}}.}

\seealso{
\code{\link{sub_model}}.
}