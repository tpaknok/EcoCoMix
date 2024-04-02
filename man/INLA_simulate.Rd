\name{INLA_simulate}
\alias{INLA_simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating responses from INLA models
}
\description{
Model diagnostics using DHARMa require simulated responses (see \code{\link[DHARMa]{createDHARMa}}). This function automatically simulate the responses from an INLA model. Currently only support models of gaussian, poisson, binomial, negative binomial (type 0), tweedie, and zero-inflated negative binomial (type 1) model.
}
\usage{
INLA_simulate(object,
              nsim = 250,
              Ntrails = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
INLA or CPR output (e.g. m$best_model)
}
  \item{nsim}{
Number of simulations
}
  \item{Ntrails}{
The number of trials for the binomial likelihood. see \code{\link[INLA]{inla}}. Defaulted as 1.
}
}
\details{
The simulated response can be passed to \code{\link[DHARMa]{createDHARMa}}, creating a DHARMa object for further diagnostic.
}
\value{
A data frame containing simulated responses, with each column representing each simulation run.
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
\code{\link[DHARMa]{createDHARMa}}
}
\examples{
library(phytools)
library(tidyverse)
library(CPR)
data("KSR")
data("KSR_MLtree")
data("KSR_EF")

VCV_sp <- vcv(KSR_MLtree) #species level phyologenetic covariance matrix using default (Brownian) model
VCV_sp <- VCV_sp[order(rownames(VCV_sp)),order(colnames(VCV_sp))]

inla_m <- CPR(formula = bugs~Real.rich+f(comm,model="generic0",Cmatrix=Phylo),
              priors=NULL,
              df = KSR_EF,
              VCV_sp = VCV_sp,
              comm=KSR,
              family="nbinomial",
              optim.lambda = T)

library(DHARMa)
sres <- INLA_simulate(inla_m$best_model)
DHARMa_res <- createDHARMa(sres,
                           KSR_EF$bugs,
                           fittedPredictedResponse=inla_m$predictedfittedresponse,
                           integerResponse = T)
plot(DHARMa_res,quantreg=F)

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
