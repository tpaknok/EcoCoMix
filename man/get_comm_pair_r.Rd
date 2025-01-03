\name{get_comm_pair_r}
\alias{get_comm_pair_r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic covariance or correlation Matrix at the community level
}
\description{
The function computes the expected covariances or correations of any community-level properties between communities.
}
\usage{
get_comm_pair_r(comm, VCV_sp, comm_kronecker = NULL, force.PD = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{
  A data frame describing species (column) abundance / presence-absencene in each site, or sample in general (row).
}
  \item{VCV_sp}{
  A phylogenetic covariance (not correlation!) matrix based on any evolutionary models.
}
  \item{comm_kronecker}{
  For internal calculation only in \code{\link[CPR]{CPR_spaMM}}.
}
  \item{force.PD}{
  Forcing Positive Definite using \code{\link[Matrix]{nearPD}}. Defaulted as FALSE
}
}
\details{
The function upscales species-level covariance to community-level covariance, which is subsequently converted to a correlation matrix.
}
\value{
  \item{corM}{Expected correlation between communities.
  }
  \item{covM}{Expected covariance between communities.
  }
  \item{comm_kronecker}{For internal uses to speed up the optimization processes.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Toby P.N. Tsang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(phytools)
library(tidyverse)
library(CPR)
library(ape)
library(DHARMa)

data(KSR)
data(KSR_MLtree)
data(KSR_EF)

phy_comm_M <- get_comm_pair_r(KSR,vcv(KSR_MLtree))$corM

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
