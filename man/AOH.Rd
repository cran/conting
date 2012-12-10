\name{AOH}
\alias{AOH}
\docType{data}
\title{
Alcohol, Obesity and Hypertension: A Complete 4 * 3 * 2 Table

}
\description{
491 subjects are cross-classified according to the three factors: hypertension (H, 2 levels), 
obesity (O, 3 levels) and alcohol (A, 4 levels). There are a total of 24 cells in the
table.
}
\usage{data(AOH)}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{Counts in each cell of table.}
    \item{\code{alc}}{A factor with levels \code{0} \code{1-2} \code{3-5} \code{6+} indicating
			the classification of alcohol intake of drinks per day.}
    \item{\code{obe}}{A factor with levels \code{low} \code{average} \code{high} indicating
			the classification of obesity.}
    \item{\code{hyp}}{A factor with levels \code{yes} \code{no} indicating the classification of
			hypertension.}
  }
}
\details{
This data are from a study in Western Australia. The study copied a larger study 
from USA. See Knuiman & Speed (1988) for more details.

For details on the function \code{\link{bcct}} applied to this data, see Overstall & King (2012).
}
\source{
Knuiman, M.W. & Speed, T.P. (1988) Incorporating Prior Information into the 
Analysis of Contingency Tables. Biometrics, \bold{44}, 1061--1071. 
}
\references{
Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\examples{
data(AOH)
summary(AOH)
}
\keyword{datasets}
