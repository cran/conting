\name{heart}
\alias{heart}
\docType{data}
\title{
Risk Factors for Coronary Heart Disease: A Complete 2^6 Table
}
\description{
1841 men are cross-classified according to six risk factors for 
coronary heart disease: smoking (A, 2 levels), strenuous mental 
work (B, 2 levels), strenuous physical work (C, 2 levels), 
systolic blood pressure (D, 2 levels), ratio of alpha and beta 
lipoproteins (E, 2 levels) and family anamnesis of coronary 
heart disease (F, 2 levels).
}
\usage{data(heart)}
\format{
  A data frame with 64 observations on the following 7 variables.
  \describe{
        \item{\code{y}}{Counts in each cell of table.}
    \item{\code{A}}{A factor with levels \code{yes} \code{no} indicating 
			smoking status.}
    \item{\code{B}}{A factor with levels \code{yes} \code{no} indicating 
			strenuous mental work.}
    \item{\code{C}}{A factor with levels \code{yes} \code{no} indicating 
			strenuous physical work.}
    \item{\code{D}}{A factor with levels \code{yes} \code{no} indicating 
			systolic blood pressure.}
    \item{\code{E}}{A factor with levels \code{yes} \code{no} indicating 
			high ratio of alpha and beta lipoproteins.}
    \item{\code{F}}{A factor with levels \code{yes} \code{no} indicating 
			a family anamnesis of coronary heart disease.}
}}
\details{
For more details on this data see Edwards & Havranek (1985).

For details on the function \code{\link{bcct}} applied to this data, see Overstall 
& King (2012).
}

\source{
Edwards, D. & Havranek, T. (1985) A fast procedure for model search 
in multidimensional contingency tables. Biometrika, \bold{72}, 339--351.
}
\references{
Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\examples{
data(heart)
summary(heart)
}
\keyword{datasets}
