\name{find_cens}
\alias{find_cens}

\title{
Find Censored Cells
}
\description{
Given all the sources and the censored source of an incomplete contingency table, this 
function will find the censored cells. 
}
\usage{
find_cens(sources, cens_source, data=NULL, unobs.level = "un", obs.level = "obs")
}

\arguments{
  \item{sources}{
An object of class \code{"formula"}, which details the sources in the incomplete contingency table.
}
  \item{cens_source}{
An object of class \code{"formula"}, which details the source which is subject to censoring in the 
incomplete contingency table.
}
  \item{data}{
An object of class \code{"data.frame"} (or \code{"table"}) containing the variables in the model. 
If the model variables are not found in \code{data}, the variables are taken from 
\code{environment(formula)}, typically the environment from which \code{find_cens} is called.
}
  \item{unobs.level}{
The character string used to label the source level corresponding to not observing the 
individuals in the cell.
}
  \item{obs.level}{
The character string used to label the source level corresponding to observing the 
individuals in the cell.
}
}
\details{
Sometimes one of the sources (termed the censored source) used to estimate closed populations 
observes individuals which are not members of the target population. In this case we assume 
that when this source observes an individual that has been observed by at least one other 
source, then it is a member of the target population. However those individuals only observed 
by the censored source contain a mixture of members of the target and non-target populations. This 
means that the observed cell count acts as an upper bound on the true cell count. For more details 
on this approach, see Overstall et al (2014) and Overstall & King (2014). This function identifies 
the cells which are censored (i.e. correspond to only being observed by the censored source).
}
\value{
The function will output a numeric vector containing the cell numbers of the censored cells. These are 
used by the \code{\link{bict}} and \code{\link{bictu}} functions.
}
\references{
Overstall, A.M., King, R., Bird, S.M., Hutchinson, S.J. & Hay, G. (2014) Incomplete contingency 
tables with censored cells with application to estimating the number of people who inject drugs 
in Scotland. \emph{Statistics in Medicine}, \bold{33 (9)}, 1564--1579.

Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Journal of Statistical Software}, \bold{58 (7)}, 
1--27. \url{http://www.jstatsoft.org/v58/i07/}
}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}

\seealso{
\code{\link{bict}}.
}
\examples{
data(ScotPWID)
## Load the ScotPWID data. In this dataset, the S4 source corresponding 
## to the HCV database is subject to censoring. We use find_cens to find
## the censored cells.
find_cens(sources=~S1+S2+S3+S4,cens_source=~S4,data=ScotPWID)
## It will produce the vector with the following elements:
##[1]   9  25  41  57  73  89 105 121
## Let's look at these cells
ScotPWID[find_cens(sources=~S1+S2+S3+S4,cens_source=~S4,data=ScotPWID),]
## It will produce:
#      y S1 S2 S3  S4 Region Gender   Age
#9   122 un un un obs    GGC   Male Young
#25  135 un un un obs    GGC   Male   Old
#41   48 un un un obs    GGC Female Young
#57   38 un un un obs    GGC Female   Old
#73  134 un un un obs   Rest   Male Young
#89  104 un un un obs   Rest   Male   Old
#105  78 un un un obs   Rest Female Young
#121  25 un un un obs   Rest Female   Old
}