\name{CPR_spaMM}
\alias{CPR_spaMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modeling community-level data controlling compositional autocorrelation.
}
\description{
The function uses fitme (package spaMM) to model compositional autocorrelation when analyzing community-level data. Phylogenetic relationship between species can be provided or optimized (assuming a Brownian motion in the latter case)
}
\usage{
CPR_spaMM(formula,
          data,
          VCV_sp,
          comm,
          optim.lambda = T,
          original.VCV = T,
          AIC_threshold = -4,
          method.spaMM = "REML",
          true_VCV = NULL,
          control.optim = NULL,
          comm_kronecker = NULL,
          int_model = T,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a formula specified in spaMM format (see \code{\link[spaMM]{fitme}}). Use "comp_id"" to represent community-level random effect controlling compositional autocorrelation (see Example).
}
  \item{data}{
  a data frame containing the variables in the model.
}
  \item{VCV_sp}{
  a covariance (not correlation!) matrix illsutrating phylogenetic relationship between species
}
  \item{comm}{
  a site(row)-species(column) matrix. Pleasure ensure species names are consistent between VCV_sp and comm.
}
  \item{optim.lambda}{
  TRUE/FALSE. Whether lambda should be optimized (through finding value with the lowest cAIC).
}
  \item{original.VCV}{
  Return model results based on the provided VCV_sp. Defaulted as TRUE.
}
  \item{AIC_threshold}{
  cAIC threshold to include compositional random effect. Defaulted as -4 (i.e., including compositional random effect needs to reduce cAIC by 4 to be considered as improving the model).
}
  \item{method.spaMM}{
  "ML" (Default), "REML", or PGL/L. Will be passed to the method argument of \code{\link[spaMM]{fitme}}.
}
  \item{true_VCV}{
  True phylogenetic relationship between species. Useful in simulations but not empirical analyses.
}
  \item{control.optim}{
  A list of control parameters for \code{\link[stats]{optim}} (See the documentation of optim when using this argument).
}
  \item{comm_kronecker}{
  Arguments used to speed up internal calculations. Defaulted as NULL.
}
  \item{int_model}{
  Obtain results for intercept-only models (both with and without the compositional random effects). Defaulted as TRUE.
}
  \item{\dots}{
  Optional arguments passed to \code{\link[spaMM]{fitme}}.
}
}
\details{
Can provide extra specifications to the fitme function to speed up calculations / obtain more reliable results. These technical details, however, are non-trivial. Please see the documentation of \code{\link[spaMM]{fitme}}, \code{\link[spaMM]{inits}}, and \link[stats]{optim}. If ignored, the model will be run using the defaults of these functions.
}
\value{
a list of 19 objects.
\item{best_model}{Results of the best model (i.e., the model with optimized lambda if its dAIC is lower than the intercept-only model by 4 (or other self-defined threshold.)}
\item{best_model_satt}{Results of the best model based on the Satterthwaite approximation, which provides better estimations of p-values.}
\item{optimized_lambda_model,optimized_lambda_model_satt}{Same as best_model and best_model_satt, but this is for the model with optimized phylogenetic structure.}
\item{original_VCV_model,original_VCV_m_satt}{Same as above, but for models based on the provided phylogenetic structure.}
\item{true_model,true_model_satt}{Same as above, but for models based on the true phylogenetic structure (only useful for simulations).}
\item{without_comp_model}{Model results without considering compositional autocorrelation.}
\item{without_comp_anova}{Model results without considering compositional autocorrelation in ANOVA format. Easier to compare with results from mixed models based on Satterthwaite approximation.}
\item{intercept_only_model}{Intercept-only model results (with compositional correlation).}
\item{star_model}{Model results without considering phylogenetic dependence between species (i.e., all species are phylogenetically independent).}
\item{AIC}{cAIC of all models above.}
\item{optimized_lambda,optimized_lambda_int}{The optimized lambda based on the full model (former) and intercept-only model(latter).}
\item{conv}{whether the following models have reached convergence: the best model, the full model based on optimized phylogenetic structure, and the model based on the provided phylogenetic relationship. Yes=1, No=0.}
}
\author{
Toby P.N. Tsang
}

\seealso{
\code{\link[spaMM]{fitme}}, \code{\link[stats]{optim}}
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

m_LAI <- CPR_spaMM(LAI~Real.rich+corrMatrix(1|comp_id),
          data=KSR_EF,
          comm=KSR,
          VCV_sp = vcv(KSR_MLtree),
          method.spaMM="REML")

m_LAI$best_model
m_LAI$best_model_satt
plot(simulateResiduals(m_LAI$optimized_lambda_model,re.form=NULL))
}
