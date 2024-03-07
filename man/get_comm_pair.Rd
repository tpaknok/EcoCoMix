\name{get_comm_pair_r}
\alias{get_comm_pair_r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Getting pairwise phylogenetic correlation between communities
}
\description{
Getting pairwise phylogenetic correlation between communities
}
\usage{
get_comm_pair_r(comm,V,force.PD = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{a species-site matrix}
  \item{V}{a phylogenetic covariance matrix at species levels}
  \item{force.PD}{find the closest positive definite matrix at community levels.
  Default as TRUE. Note that if matrix is not positive definite subsequent analyses can fail.}
}

\details{
  The resulting phylogenetic correlation matrix is based on the additive law of covariance.
  Imagine two communities with species (X,Y) and species (X,Z). Covariance at community levels   is Cov(X+Y,X+Z) = Var(X)+Cov(X,Z)+Cov(X,Y)+Cov(Y,Z). All of these can be obtained from the
  species level phylogenetic matrix.
}
\value{
  a phylogenetic correlation matrix at community levels.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Toby P.N. Tsang
}

\examples{
library(phytools)

}
