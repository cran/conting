\name{ScotPWID}
\alias{ScotPWID}
\docType{data}
\title{
People Who Inject Drugs in Scotland 2006: An Incomplete 2^7 Table
}
\description{
5670 people who inject drugs (PWID) in Scotland in 2006 are osberved by 
four sources: social enquiry reports (\code{S1}), hospital records (\code{S2}), 
Scottish drug misuse database (\code{S3}) and Hepatitis C virus  (HCV) diagnosis 
database (\code{S4}). The PWID are further cross-classified according to three 
additional factors: region (\code{Region}; 2 levels), gender (\code{Gender}; 2 
levels) and age (\code{Age}; 2 levels).
}
\usage{data(ScotPWID)}
\format{
  A \code{"data.frame"} with 128 observations on the following 8 variables.
  \describe{
    \item{\code{y}}{
Counts in each cell of the table with NAs for the cells corresponding to not being observed by any of the sources.
}
    \item{\code{S1}}{
A factor with levels \code{un} \code{obs} indicating whether source S1 observed the PWID.
}
    \item{\code{S2}}{
A factor with levels \code{un} \code{obs} indicating whether source S2 observed the PWID.
}
    \item{\code{S3}}{
A factor with levels \code{un} \code{obs} indicating whether source S3 observed the PWID.
}
    \item{\code{S4}}{
A factor with levels \code{un} \code{obs} indicating whether source S4 observed the PWID.
}
    \item{\code{Region}}{
A factor with levels \code{GGC} \code{Rest} indicating 	the region (\code{GGC} = Greater Glasgow & Clyde, \code{Rest} = Rest of Scotland).
}
    \item{\code{Gender}}{
A factor with levels \code{Male} \code{Female} indicating gender.
}
    \item{\code{Age}}{
A factor with levels \code{Young} \code{Old} indicating age (\code{Young} = <35 years, \code{Old}=35+ years).
}
  }
}
\details{
Note that the PWID observed by source \code{S4}, the HCV database, are not necessarily current 
PWID. They are people who have a history of drug use. Therefore the count in the cell 
corresponding to only being observed by the HCV database is an overcount. Overstall et al 
(2014) use a modelling approach whereby the count in the cell corresponding to only being 
observed by the HCV database is missing and the observed value acts as an upper bound. For 
more details on the dataset see King et al (2013).

For details on the function \code{\link{bict}} applied to this data, see Overstall & King (2014).
}
\source{
King, R., Bird, S. M., Overstall, A. M., Hay, G. & Hutchinson, S. J. (2013) Injecting drug
users in Scotland, 2006: Listing, number, demography, and opiate-related death-rates.
\emph{Addiction Research and Theory}, \bold{21}, 235-246.
}
\references{
Overstall, A.M., King, R., Bird, S.M., Hutchinson, S.J. & Hay, G. (2014) Incomplete contingency 
tables with censored cells with application to estimating the number of people who inject drugs 
in Scotland. \emph{Statistics in Medicine}, \bold{33}, 1564--1579.

Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\examples{
data(ScotPWID)
summary(ScotPWID)
}
\keyword{datasets}
