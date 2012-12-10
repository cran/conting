\name{rabbit}
\alias{rabbit}
\docType{data}
\title{
Cottontail Rabbits: An Incomplete 2^3 Table
}
\description{
115 cottontail rabbits are observed by three different sources: first live 
trapping (first), second live trapping (second) and hunter's bag (hunter).
}
\usage{data(rabbit)}
\format{
  A data frame with 8 observations on the following 4 variables.
  \describe{
    \item{\code{first}}{
A factor with levels \code{obs} \code{un} indicating whether this source 
			observed the rabbit.}
    \item{\code{second}}{
A factor with levels \code{obs} \code{un} indicating whether this source 
			observed the rabbit.}
    \item{\code{hunter}}{
A factor with levels \code{obs} \code{un} indicating whether this source 
			observed the rabbit.}
    \item{\code{y}}{
Counts in each cell of table with an NA for the cell corresponding to not 
			being observed by any of the sources.}
  }
}
\details{
Between 1965 and 1970, the Rose Lake Wildlife Research Station in Michigan, USA collected counts of 
cottontail rabbits so they could measure the population. In each year, the rabbits are trapped for 
about three weeks before the hunting season, marked and then returned to the trapping area. There are 
two trapping periods, denoted first and second, and these are the first two sources. The final source 
is the hunter's bag in the hunting season. The data here is from 1967. For more details and the full 
dataset see Bishop et al (1975, pg 245).

For details on the function \code{\link{bict}} applied to this data, see Overstall & King (2012).
}
\source{
Bishop, Y.M.M., Fienberg, S.E. & Holland, P.W. (1975) \emph{Discrete Multivariate Analysis}. MIT Press.
}
\references{
Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\examples{
data(rabbit)
summary(rabbit)
}
\keyword{datasets}
