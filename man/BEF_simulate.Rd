\name{BEF_simulate}
\alias{BEF_simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating a BEF relationship and evaluting different analytical methods}
\description{
Simulating a BEF relationship with/without phylogenetic signals, followed by evaluating different analytical methods based on type-I/type-II errors, estimated phylogenetic signals, and R2.
}
\usage{
BEF_simulate(comm,
             CV_sp = NULL,
             nspp=15,
             nsite=50,
             min_richness=1, max_richness=4,
             spaMM_formula,
             b0=0,
             b1=0,
             signals_X="phy_cor",
             noise_mean = 0,
             noise_sd = 1,
             non_phy_cor_mean = 0,
             non_phy_cor_sd = 1,
             lambda_true=1,
             scale_all=F,
             conv_fail_drop = T,
             ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{A species-site matrix.}
  \item{VCV_sp}{A phylogenetic covariance matrix at species levels. If ignored, a stochastic birth-death tree will be generated.}
  \item{nspp}{Species pool size. The default is 15.}
  \item{nsite}{Number of sites. The default is 50.}
  \item{min_richness}{The minimum species richness among all sites. The default is 1.}
  \item{max_richness}{The maximum species richness among all sites. The default is 1.}
  \item{spaMM_formula}{A formula in spaMM syntax. See Details for specifying phylogenetic random effects}
  \item{b0}{The value of b0, which is the global intercept. The default is 0}
  \item{b1}{The value of b1, which is the effect of X on Y.}
  \item{signals_X}{Type of X simulated or used. If signals_X = "sr", species richness of simulated communities will be used. If signals_X = "phy_cor", a continuous variable with phylogenetic signals will be generated. If signals_X = "no_phy_cor", a continuous variable without phylogenetic signals will be generated. The mean and S.D. of the normal distribution should be specified in argument non_phy_cor_mean and non_phy_cor_sd. The default is "phy_cor".}
  \item{noise_mean}{The mean of the random noise generated as a normal distribution. The default is 0}
  \item{noise_sd}{The S.D. of the random noise generated as a normal distribution. The default is 1}
  \item{non_phy_cor_mean}{The mean of the normal distribution used to generate X if it is specified as "no_phy_cor". The default is 0.}
  \item{non_phy_cor_sd}{The S.D. of the normal distribution used to generate X if it is specified as "no_phy_cor". The default is 1.}
  \item{lambda_true}{Phylogenetic signals in how species contribution to ecosyste functions. This can lead to phylogeny influencing compositional autocorrelation between communities. The default is 1.}
  \item{scale_all}{Standardizing all variables after centering, The default is FALSE.}
  \item{conv_fail_drop}{Drop results when spaMM failed to converge. The default is TRUE.}
  \item{...}{Othert arguments that will be passed to CPR_spaMM.}
}
\details{
The simulated relationship is \deqn{Y = X_{0}+b_{1}X_{1}+X_{2}+X_{3}}

\eqn{X_{0}} is the population-level intercept. \eqn{X_{2}} is a variable with compoisitional autocorrelation, with the extent of phylogenetic signal depending on lambda_true. \eqn{X_{1}} can be species richness, a variable with or without compositional autocorrelation, depending on signals_X. \eqn{b_{1}} is the poluation-level effect of \eqn{X_{1}} and can be modified by users. \eqn{X_{3}} is the random noise, defaulted as \eqn{N(0,1)}.

\eqn{X_{2}} (and \eqn{X_{1}} if signals_X = "phy_cor") follows \deqn{MVN(0,{C_{comm})}}

\eqn{C_{comm}} represents the correlation matrix between communities. The estent of phylogenetic influence depends on lambda_true.
}

\value{
A list of data frames with the length equals the number of simulation.

The following indicates whether X1 is statistically significant based on different model specification. 1=Y; 0=N.
\item{m_optim_sig}{Model based on the optimized phylogenetic structure.}
\item{m_true_sig}{Model based on the phylogenetic structure.}
\item{m_original_sig}{Model based on the phylogenetic structure.}
\item{m_best_sig}{The most parsimonious model.}
\item{m_without_comp_sig}{Model without considering the compositional autocorrelation.}

The following indicates the estimated slope of X1:
\item{m_optim_slope}{Model based on the optimized phylogenetic structure.}
\item{m_true_slope}{Model based on the true phylogenetic structure.}
\item{m_original_slope}{Model based on the provided phylogenetic structure.}
\item{m_best_slope}{The most parsimonious model.}
\item{m_without_comp_slope}{Model without considering the compositional autocorrelation.}

Optimized lambda:
\item{optim_lambda}{Optimized lambda for the model m_optim}
\item{optim_lambda_int}{Optimized lambda for the intercept-only model.}

The AIC of models:
\item{AIC_without_comp}{Model without considering compositional autocorrelation.}
\item{AIC_original_VCV}{Model based on the provided phylogenetic structure}
\item{AIC_optim}{Model based on the optimized phylogenetic structure.}
\item{AIC_star}{Model based on the star phylogeny.}
\item{AIC_true}{Model based on the true phylogenetic structure}
\item{AIC_optim_int}{Intercept-only model based on the optimized phylogenetic structure}

Simulation properties:
\item{min_richness}{The minimum richness of communities for the simulation.}
\item{max_richness}{The maximum richness of communities for the simulation.}
\item{nspp}{Species pool size.}
\item{true_lambda}{True lambda of the phylogenetic structure}
\item{r_x1x2}{Correlation between \eqn{X_{1}} and \eqn{X_{2}}}
\item{b1}{True value of \eqn{b_{1}}.}
\item{Count}{Number of iterations. >1 if any model failed to converge and if conv_fail_drop = T.}
\item{signals_X}{Input for the signals_X argument.}
\item{phylo_structure}{"Simulated" (if VCV_sp is null) or "Provided" (if VCV_sp is not null).}

Other results:
\item{optim_r2m}{Marginal R-squared of the optimized model.}
\item{optim_r2c}{Conditional R-squared of the optimized model.}
\item{NumDF}{Numerator degree-of-freedom based on the Satterthwaite approximation for the optimized model.}
\item{DenDF}{Conditional degree-of-freedom based on the Satterthwaite approximation for the optimized model.}


}

\author{
Toby P.N. Tsang
}
\seealso{
\code{\link[CPR]{CPR_spaMM}}
}
\examples{

library(phytools)
library(tidyverse)
library(CPR)
library(ape)

nspp <- 14
sim <- 50
b1 <- 0 #assume no effects of species richness
signals_X <- "sr"
lambda_true <- runif(sim) #random phylogenetic signals

count <- 1

result <- list()

spaMM_formula <-  y~x1+corrMatrix(1|comp_id)

for (i in 1:sim) {
  message("sim_",i)

  result[[i]] <- tryCatch(BEF_simulate(comm = NULL,
                                           VCV_sp = NULL,
                                           nspp=nspp,
                                           nsite=88,
                                           min_richness=1,
                                           max_richness= 4,
                                           spaMM_formula=spaMM_formula,
                                           b1=b1,
                                           signals_X=signals_X,
                                           noise_mean = 0,
                                           noise_sd = 0.01, #very low noise
                                           lambda_true= lambda_true[[i]],
                                           conv_fail_drop = T,
                                           scale_all=F,
                                           control.optim=list(factr=1e12),
                                           optim.lambda=T,
                                           init=list(), #specification for spaMM - calculations will be faster
                                           int_model=F)
                              error = function(e) e)
  }


all_result <- as.data.frame(do.call(rbind,result))
}
