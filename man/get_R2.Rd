\name{get_R2}
\alias{get_R2}

\title{
Provide marginal and conditional R-squared for outputs from spaMM
}
\description{
Provide marginal and conditional R-squared for outputs from spaMM
}
\usage{
get_R2(spaMM_m,spaMM_m_int_only=NULL)
}

\arguments{
  \item{spaMM_m}{model output from spaMM
  }
  \item{spaMM_m_int_only}{intercept-only model from spaMM. Only needed for negative binomial distribution.
  }
}

\details{
The function is based on \code{\link[performance]{R2}} from the performance package, which does not support spaMM. Note that the model currently only supports gassuian and negative binomial (negbin2) distribution. If the model doesn't contain any compositional autocorrelation, it's probably easier to use other packages or functions (e.g. lme4, glmmTMB, lm, glm), which is better supported by other R packages.

You can also calculate the correlation coefficient between observed and predicted responses (prediction can be based on fixed effects only or full model, depending on the purposes). This is technically not the same as marginal and conditional R-squared, but should be a easy solution if you have models based on distributions not supported here.
}
\value{
  \item{R2m}{Marginal R-squared.
  }
  \item{R2c}{Conditional R-squared.
  }
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
library(EcoCoMix)
library(ape)
library(DHARMa)

data(KSR)
data(KSR_MLtree)
data(KSR_EF)

### Normal distribution
m_LAI <- EcoCoMix(LAI~Real.rich+corrMatrix(1|comp_id),
                   data=KSR_EF,
                   comm=KSR,
                   VCV_sp=vcv(KSR_MLtree),
                   method.spaMM="REML")

get_R2(m_LAI$best_model) #get R2 for best model
get_R2(m_LAI$original_VCV_model) #get R2 for model based on the provided phylogenetic covariance matrix between species

### nbinom2
m_bugs <- EcoCoMix(bugs~Real.rich+corrMatrix(1|comp_id),
                   data=KSR_EF,
                   comm=KSR,
                   VCV_sp = vcv(KSR_MLtree),
                   family="negbin2")

# for the full model with lambda = 1
# obtain phylogenetic correlation matrix at the community level
phy_comm_M <- get_comm_pair_r(KSR,vcv(KSR_MLtree))$corM

m_bugs_int <- fitme(bugs~1+corrMatrix(1|comp_id),corrMatrix=as_precision(phy_comm_M),data=KSR_EF,family="negbin2")

get_R2(m_bugs$original_VCV_model,m_bugs_int)

# for the full model with optimized lambda (also the best model in this case)

# modify the species level phylogenetic covariance matrix based on the optimized
optim_VCV_sp <- vcv(KSR_MLtree)*m_LAI$optimized_lambda  lambda
diag(optim_VCV_sp) <- diag(vcv(KSR_MLtree))

# obtain phylogenetic correlation matrix at the community level
phy_comm_M <- get_comm_pair_r(KSR,optim_VCV_sp)$corM

m_bugs_int <- fitme(bugs~1+corrMatrix(1|comp_id),corrMatrix=as_precision(phy_comm_M),data=KSR_EF,family="negbin2")

get_R2(m_bugs$best_model,m_bugs_int)

### best model without random effect
m_litter <- EcoCoMix(litter2012~Real.rich+corrMatrix(1|comp_id),
                   data=KSR_EF,
                   comm=KSR,
                   VCV_sp = vcv(KSR_MLtree))

m_litter_best <- fitme(litter2012~Real.rich,data=KSR_EF) #fitme approach
pseudoR2(m_litter_best)

m_litter_best2 <- lm(litter2012~Real.rich,data=KSR_EF)
summary(m_litter_best2) #same result
}
