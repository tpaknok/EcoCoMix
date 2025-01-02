\name{get_R2}
\alias{get_R2}

\title{
Provide marginal and conditional R-squared for outputs from spaMM
}
\description{
Provide marginal and conditional R-squared for outputs from spaMM
}
\usage{
get_R2(model)
}

\arguments{
  \item{model}{model output from spaMM
}
}
\details{
The function is based on \code{performance}{r2}, which does not support spaMM. Note that the model currently only supports gassuian distribution. If the model doesn't contain any compositional autocorrelation, it's probably easier to use other packages (e.g. lme4, glmmTMB), which is better supported by other R packages and the calculation of R-squared is easier.

Another solution is to calculate the correlation coefficient between observed and predicted responses (prediction can be based on fixed effects only or full model). This is technically not the same as marginal and conditional R-squared, but should be a easy solution.

}
\value{
\item{R2m}{Marginal R-squared.}
\item{R2c}{Conditional R-squared.}
}

\author{
Toby P.N. Tsang
}

\seealso{
\code{performance}{r2}
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
get_R2(m_LAI$best_model)
