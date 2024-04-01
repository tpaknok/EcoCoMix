\name{BEF_simulate_eval}
\alias{BEF_simulate_eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate simulation results
}
\description{
Evaluating GLS (CPR_GLS) and INLA (CPR) results analyzing simulated data from BEF_simulate. The function will evaluate type-I / type-II errors, slope estimates, AIC/wAIC, and the optimized lambda. Results of each simulation will also be provided.
}
\usage{
BEF_simulate_eval(models, type = "INLA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
a list of outputs from CPR_GLS or CPR
}
  \item{type}{
"INLA" or "GLS". Defaulted as "INLA"
}
}
\details{
If type = "GLS", AIC will be calculated. If type = "INLA", wAIC will be calculated.

The function always evaluate the proportion of simulate runs producing a statistically significant relationship (i.e., p < 0.05). Thus whether this represents type-I error / power (1-type-II error) will depend on the true b1 specified in BEF_simulate.

}
\value{
  \item{summary_stat}{
  return the proportion of significant results, min max and mean slope estimates across runs.
  }
  \item{result}{
  detailed results of each simulation run on the slope estimates of X.
  }
  \item{AIC}{
  AIC (GLS) or wAIC (INLA) of each model. The function evaluates the model with optimized lambda, the model without considering phylogeny, and the model with the user-supplied phylogenetic covariance matrix.
  }
  \item{optim_lambda}{
  Sumamry statistics of optimized lambda.
  }
  \item{type}{
  INLA/GLS
  }
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
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

VCV_sp <- vcv(KSR_MLtree)

## simulation setup
sim <- 20  #20 runs in total
optim_sim <- 10 #optimizing lambda in ten runs only to minimize computational time
optim_seq <- sample(c(rep(TRUE, optim_sim), rep(FALSE, sim-optim_sim)), sim ,replace = F)

sim_data_lambda_1 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0,
                                  signals_X = "sr",
                                  signals_intercept = TRUE,
                                  signals_slope = FALSE,
                                  lambda_true = 1,
                                  sim = sim)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
GLS_formula <- y~x
##

GLS_lambda_1 <- lapply(1:sim, function(x) CPR_GLS(formula = GLS_formula,
                                                  df = sim_data_lambda_1$sim_dat[[x]],
                                                  VCV_sp = VCV_sp,
                                                  optim.lambda = optim_seq[[x]],
                                                  comm = KSR))

INLA_lambda_1 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                             priors=NULL,
                                             df = sim_data_lambda_1$sim_dat[[x]],
                                             VCV_sp = VCV_sp,
                                             comm=KSR,
                                             family= "gaussian",
                                             inla.rerun=3,
                                             optim.lambda = optim_seq[[x]]))

GLS_eval_lambda_1 <- BEF_simulate_eval(GLS_lambda_1,type="GLS")
INLA_eval_lambda_1 <- BEF_simulate_eval(INLA_lambda_1,type="INLA")
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
