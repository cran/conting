\name{spina}
\alias{spina}
\docType{data}
\title{
Persons born with Spina Bifida: An Incomplete 2 * 2 * 2 * 3 Table 
}
\description{
621 people born with Spina Bifida (a congenital disorder) in the state 
of New York between 1969 and 1974 are observed by three sources: 
birth certificates (\code{S1}); death certificates (\code{S2}); and medical 
rehabilitation lists (\code{S3}). The people are also cross-classified 
according to their ethnicity (\code{eth}; 3 levels). 
}
\usage{data(spina)}
\format{
  A \code{"data.frame"} with 24 observations on the following 5 variables.
  \describe{
    \item{\code{y}}{
Counts in each cell of the table with NAs for the cells corresponding to not being observed by any of the sources.
}
    \item{\code{S1}}{
A factor with levels \code{un} \code{obs} indicating whether the birth certificate source observed the person.
}
    \item{\code{S2}}{
A factor with levels \code{un} \code{obs} indicating whether the death certificate source observed the person.
}
    \item{\code{S3}}{
A factor with levels \code{un} \code{obs} indicating whether the medical rehabilitation source observed the person.
}
    \item{\code{eth}}{
A factor with levels \code{afro-american} \code{caucasian} code{other} indicating the ethnicity of the person (\code{afro-american} = Afro-American, 
\code{causcasian} = Caucasian, \code{other} = Other).
}
  }
}
\details{
See Madigan & York (1997), and the references therin, for more details on the study.

For details on the function \code{\link{bict}} applied to this data, see Overstall & King (2014).
}
\source{
Madigan, D. & York, J.C. (1997) Methods for Estimation of the Size of a Closed Population. 
\emph{Biometrika}, \bold{84}, 19--31.
}
\references{
Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\examples{
data(spina)
summary(spina)
}
\keyword{datasets}
