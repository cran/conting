\name{summary.bcct}
\alias{summary.bcct}
\alias{print.sbcct}
\alias{summary.bict}
\alias{print.sbict}

\title{
Summary of \code{bcct} and \code{bict} Objects
}

\description{
These functions produce summaries of objects of class \code{"bcct"} and \code{"bict"}. They also control 
how these summaries are printed.
}
\usage{
\method{summary}{bcct}(object, n.burnin = 0, thin = 1, cutoff = 0.75, statistic = "X2", 
best = NULL, scale = 0.1, prob.level = 0.95, ...)

\method{print}{sbcct}(x, ..., digits = max(3, getOption("digits") - 3))

\method{summary}{bict}(object, n.burnin = 0, thin = 1, cutoff = 0.75, statistic = "X2", 
best = NULL, scale = 0.1, prob.level = 0.95, ...)

\method{print}{sbict}(x, ..., digits = max(3, getOption("digits") - 3))

}

\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"}.
}
  \item{x}{
An object of class \code{"sbcct"} or \code{"sbict"} produced as a result of a call to the functions \code{summary.bcct} or 
\code{summary.bict}, respectively.
}
  \item{n.burnin}{
An optional argument giving the number of iterations to use as burn-in. 
The default value is 0.
}
  \item{thin}{
An optional argument giving the amount of thinning to use, i.e. the computations are 
based on every \code{thin}-th value in the MCMC sample. The default value is 1, i.e. no
thinning.
}
  \item{cutoff}{
An optional argument giving the cutoff posterior probability for displaying posterior 
summary statistics of the log-linear parameters. Only those log-linear parameters with 
a posterior probability greater than \code{cutoff} will be returned as part of the output. 
The default value is 0.75.
}
  \item{statistic}{
An optional argument giving the discrepancy statistic to use for calculating the Bayesian p-value. It can be one of 
\code{c("X2","FreemanTukey","deviance")} which correspond to the different statistics: 
\code{"X2"} = Chi-squared statistic, \code{"FreemanTukey"} = Freeman-Tukey statistic, 
\code{"deviance"} = deviance statistic. See Overstall & King (2014), and references 
therein, for descriptions of these statistics.
}
  \item{best}{
An optional argument for controlling how the posterior model probabilities are returned 
as output. The function will return details on the \code{best} models with the highest posterior 
model probabilities. The default value is \code{NULL}. If not \code{NULL} than this argument takes precedent 
over \code{scale}.
}
  \item{scale}{
An optional argument for controlling how the posterior model probabilities are returned 
as output. The function will return details on the models with the posterior model probability 
larger than \code{scale} times the probability of the posterior modal model. The default value is
0.1.
}
  \item{prob.level}{
An optional argument giving the probability content of the highest posterior density intervals (HPDIs). 
The default value is 0.95.
}
  \item{digits}{
An optional argument controling the rounding of output.
}  
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\details{The functions \code{summary.bcct} and \code{summary.bict} rely on the functions 
\code{\link{inter_stats}}, \code{\link{mod_probs}}, \code{\link{bayespval}}, and (in the case 
of \code{summary.bict}) \code{\link{total_pop}}. For extra information about the output from these 
functions, see the associated help files.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
comuter memory storage.
}

\value{
The function \code{summary.bcct} will return an object of class \code{"sbcct"} which is a list 
with the following components.

\item{BETA}{An \code{n.sample} by p matrix containing the sampled values of the log-linear parameters, where p 
is the number of log-linear parameters in the maximal model. For elements of this matrix which correspond 
to a log-linear parameter which is not present for the current model a zero is returned.}

\item{MODEL}{A vector of length \code{n.sample} giving the samlpled model indicators in hexadecimal format.}

\item{SIG}{A vector of length \code{n.sample} giving the sampled values for sigma^2 under the Sabanes-Bove & Held prior. If the 
unit information prior is used then the components of this vector will be one.}

\item{rj_acc}{A binary vector of the same length as the number of reversible jump moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{mh_acc}{A binary vector of the same length as the number of Metropolis-Hastings moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{priornum}{A numeric scalar indicating which prior was used: 1 = \code{"UIP"}, 2 = \code{"SBH"}.}

\item{maximal.mod}{An object of class \code{"glm"} giving the fit of the maximal model.}

\item{IP}{A p by p matrix giving the inverse of the prior scale matrix for the maximal model.}

\item{eta.hat}{A vector of length n (number of cells) giving the posterior mode of the linear predictor under 
the maximal model.}

\item{save}{The argument \code{save}.}

\item{name}{The argument \code{name}.}

\item{int_stats}{A list which contains the same components as an object of class \code{"interstat"}, i.e. 
summary statistics for the log-linear parameters, see \code{\link{inter_stats}}.}

\item{mod_stats}{A list which contains the same components as an object of class \code{"modprobs"}, i.e. 
summary statistics for the posterior model probabilities, see \code{\link{mod_probs}}.}

\item{pval_stats}{A list which contains the same components as an object of class \code{"pval"}, i.e. 
summary statistics for the posterior model probabilities, see \code{\link{bayespval}}.}

The function \code{summary.bict} will return an object of class \code{"sbict"} which is a list 
with the same components as an object of class \code{"sbcct"} and the following additional 
components.

\item{Y0}{An \code{n.sample} by k matrix giving the sampled values of the missing and censored cell counts, where k is 
the total number of missing and censored cell counts.}

\item{tpop_stats}{A list which contains the same components as an object of class \code{"totpop"}, i.e. 
posterior summary statistics for the total population, see \code{\link{total_pop}}.}

The functions \code{\link{print.sbcct}} and \code{\link{print.sbict}} will print out the MCMC acceptance rates, posterior 
summary statistics for the log-linear parameters, the posterior model probabilities, the Bayesian p-value 
and (in the case of \code{\link{print.sbict}}) posterior summary statistics for the total population size.

}

\references{
Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Journal of Statistical Software}, \bold{58 (7)}, 
1--27. \url{http://www.jstatsoft.org/v58/i07/}
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\note{
For examples see the help files for \code{\link{bcct}} and \code{\link{bict}}.
}

\seealso{
\code{\link{bcct}},
\code{\link{bict}},
\code{\link{accept_rate}},
\code{\link{bayespval}},
\code{\link{inter_stats}},
\code{\link{mod_probs}},
\code{\link{total_pop}}.

}
