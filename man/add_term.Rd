\name{add_term}
\alias{add_term}
\alias{drop_term}
\alias{prop_mod}
\title{
Determines Model Moves Given Current Model
}
\description{
These functions are used to detemine which models we can propose moves to, given the 
current model in the MCMC algorithm, and the principle of marginality.
}
\usage{
add_term(curr.index, data, maximal.mod)
drop_term(curr.index, data, maximal.mod)
prop_mod(curr.index,data,maximal.mod,null.move.prob=0.5)

}

\arguments{
  \item{curr.index}{
A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the current model.
}
  \item{data}{
An object of class \code{"data.frame"} containing the variables in the model.
}
  \item{maximal.mod}{
An object of class \code{"glm"} giving the fit of the maximal model.
}
  \item{null.move.prob}{
An optional scalar argument giving the probability of performing a null move, i.e. proposing a move to the 
current model. The default value is 0.5.
}
}
\details{
In the reversible jump algorithm we propose a move to a model given the current model. The function \code{prop_mod} 
implements a scheme whereby only local proposals are made, i.e. either a term is added or dropped. These 
types of move are called birth and death moves, respectively, by Forster et al (2012).

When a term is either added or dropped, we preserve the principle of marginality, e.g. we can 
only propose to add a three-way interaction if all the possible two-way interactions between the three factors 
are included in the present model.

The functions \code{add_term} and \code{drop_term} determine which terms can be added or dropped whilst preserving the 
principle of marginality.

The function \code{prop_mod} will call \code{add_term} and \code{drop_term} thus determining which terms can be added or 
dropped. With probability \code{null.move.prob} it will choose to remain in the current model; otherwise it 
will choose one of the possible terms to add or drop.
}
\value{
The functions \code{add_term} and \code{drop_term} will output a character 
vector containing the names of terms that can be dropped.

The function \code{prop_mod} will return a list with the following components.

\item{new.index}{A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the proposed model.}
\item{type}{A character string which will be one of \code{c("null","drop","add")} 
depending on the type of move proposed.}
\item{total.choices}{If \code{type} is not \code{"null"}, then \code{total.choices} will be scalar giving the total number of non-null moves available. If \code{type} is equal to \code{"null"}, then \code{total.choices} will be 0.}
\item{null.move.prob}{A scalar giving the probability of a null move.}
}
\references{
Forster, J.J., Gill, R.C. & Overstall, A.M. (2012) Reversible jump methods for generalised linear models and generalised 
linear mixed models. \emph{Statistics and Computing}, \bold{22 (1)}, 107--120.
}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\note{
This function will not typically be called by the user.
}

\examples{

data(AOH)
## Load the AOH data

maximal.mod<-glm(formula=y~(alc+hyp+obe)^3,data=AOH,x=TRUE,y=TRUE,
contrasts=list(alc="contr.sum",hyp="contr.sum",obe="contr.sum"))
## Set up the maximal model which in this case is the saturated model.

curr.index<-formula2index(big.X=maximal.mod$x,formula=y~alc+hyp+obe+hyp:obe,data=AOH)
## Set up the binary vector for the model containing all main effects and the 
## hyp:obe interaction.

add_term(curr.index=curr.index,data=AOH,maximal.mod=maximal.mod)
## See what terms we can add - will get: 

#[1] "alc:hyp" "alc:obe"

drop_term(curr.index=curr.index,data=AOH,maximal.mod=maximal.mod)
## See what terms we can drop - will get: 

#[1] "hyp:obe"

set.seed(4)
## Set the seed for reproducibility.

prop_mod(curr.index=curr.index,data=AOH,maximal.mod=maximal.mod)
## Propose a model. Will be a drop move, proposing the independence model by 
## dropping the hyp:obe interaction. The total.choices object is 3, i.e. one 
## drop move and two add moves. Specifically:

#$new.index
# [1] 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#$type
#[1] "drop"
#
#$total.choices
#[1] 3
#
#$null.move.prob
#[1] 0.5

}
