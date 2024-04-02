\name{BEF_simulate}
\alias{BEF_simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating a BEF relationship}
\description{
Simulating a BEF relationship with/without phylogenetic signals. Broadly, this can also be used to simulate any relationship with phyloogenetic signals.
}
\usage{
BEF_simulate(comm,
V,
ef_mean,
sd,
b1 = 0,
signals_X = "sr",
signals_intercept = TRUE,
signals_slope = FALSE,
lambda_true = 1,
sim = 500,
seed=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{A species-site matrix.}
  \item{V}{A phylogenetic covariance matrix at species levels.}
  \item{ef_mean}{The mean value of the simulated per-capita effect on Y of each species. The default is zero. Used only when signals_X = "phy_cor" (see Details)}
  \item{sd}{The standard deviation of the simulated per-capita effect on Y of each species. The default is one. Used only when signals_X = "phy_cor" (see Details)}
  \item{b1}{The value of b1, which is the effect of X on Y.}
  \item{signals_X}{Type of X simulated or used. If signals_X = "sr", species richness will be used. If signals_X = "phy_cor", a continuous variable with phylogenetic signals will be generated. If signals_X = "no_phy_cor", a continuous variable without phylogenetic signals (\eqn{N(0,1)}) will be generated. The default is "sr".}
  \item{signals_intercept}{Phylogenetic signals in the intercept of each community. The default is TRUE.}
  \item{signals_slope}{Phylogenetic signals in the slope of each community. The default is FALSE.}
  \item{lambda_true}{Lambda used in simulating the data. The default is 1. Ignored if signals_intecept and signals_slope are FALSE.}
  \item{sim}{Number of simulations. The default is 500.}
  \item{seed}{The seed used to simulate the data. The default is 1000.}
}
\details{
The simulated relationship is \eqn{Y = (b_{0}+Z_{1,comm})+(b_{1}+Z_{0,comm})X},

where \eqn{b_{0}} is the population-level intercept, \eqn{b_{1}} is the poluation-level effect of \eqn{X}, \eqn{Z_{0,comm}} is the random intercept for each community, \eqn{Z_{1,comm}} is the random slope for each community. If signals_intercept and signals_slope are FALSE, both random effects will be zero.

If signals_X = "phy_cor", X will be \eqn{chol({C_{comm})U}}, where \eqn{U} is a data vector randomly sampled from a Normal distirbution \eqn{N(ef\_ mean,sd)}, and \eqn{chol({C_{comm})}} is the upper triangular of the Cholesky decomposition of the phylogenetic covariance matrix at community levels. This is also the same way how random intercepts and slopes are generated
}
\value{
A list of data frames with the length equals the number of simulation. The data frames contained two columns X and Y.
}
\author{
Toby P.N. Tsang
}
\seealso{
\code{\link[CPR]{BEF_simulate_eval}}
}
\examples{
library(phytools)
library(tidyverse)
library(CPR)
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

VCV_sp <- vcv(KSR_MLtree)

sim_data_lambda_1 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0,
                                  signals_X = "sr",
                                  signals_intercept = TRUE,
                                  signals_slope = FALSE,
                                  lambda_true = 1,
                                  sim = 500)
}
