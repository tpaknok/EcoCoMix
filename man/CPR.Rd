\name{CPR}
\alias{CPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetic bayesian generalized linear mixed models for community data
}
\description{
Conducting Bayesian generalized linear mixed models for community data using INLA.
}
\usage{
CPR(formula,
    priors = NULL,
    df,
    VCV_sp,
    comm,
    family,
    optim.lambda = TRUE,
    original.VCV = TRUE,
    optim.control = NULL,
    control.compute = list(waic = TRUE, config = TRUE),
    prediction.df = "auto",
    predictedfittedresponse = "best_m",
    inla.rerun = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an INLA formula.
}
  \item{priors}{
prior specifications. If NULL default priors will be used. Use \code{inla.list.models("prior")} to see all prior available.
}
  \item{df}{
a data frame with the response and predictors
}
  \item{VCV_sp}{
a phylogenetic covariance matrix at species levels
}
  \item{comm}{
a species site matrix
}
  \item{family}{
error distribution in the glmm. Defaulted as "gaussian". Use \code{inla.list.models("likelihood")} to see all distribution available.
}
  \item{optim.lambda}{
TRUE/FALSE. If TRUE, the function will optimize the (\eqn{\lambda}).
}
  \item{original.VCV}{
TRUE/FALSE. If TRUE, the function will build a model based on the phylogenetic covariance matrix supplied by the user.
}
  \item{optim.control}{
additional arguments passing to \code{\link[stats]{optim}}.
}
  \item{control.compute}{
arguments for control.compute in the \code{\link[INLA]{inla}} function.
}
  \item{prediction.df}{
model predictions. If "auto", the model will provide predictions of a focal relatonship (e.g. ecosystem function and biodiversity) by holding other covariates at their mean (continuous) or mode (nominal) values. Users can also provide a data frame and ask INLA to make model predictions. The data frame must have the same structure as the analyzed datset, and users can define the conditions of different X while setting the response as NA.
}
  \item{predictedfittedresponse}{
compute the median posterior predictions of the best model ("best_m"), models with the phylogenetic covariance matrix supplied by users ("original_VCV"), or models without phylogenetic effect ("no_phylo"). The results can then passed to \code{\link[DHARMa]{createDHARMa}}. Defaulted as "best_m".
}
  \item{inla.rerun}{
number of times to rerun INLA models to ensure stable results. Defaulted as 1.
}
  \item{\dots}{
additional arguments passing to \code{\link[INLA]{inla}}
}
}
\details{
When specifying the random effect, use \code{Cmatrix = Phylo} to indicate the phylogenetic covariance between communities, and \code{comm} as the identity of each community. See examples. You don't need to create them in the data frame. CPR will automatically create these for you in the internal environment.

The model will always compare the wAIC between the model without considering phylogeny and the model with optimized (\eqn{\lambda}). The model without phylogeny will be selected as "best model" if the wAIC does not reduce more than 2 after including the optimized/original \eqn{\lambda}, which means \deqn{wAIC_{optimized/original} - wAIC_{without\_phylogeny} > -2}

Note that even if optimized (\eqn{\lambda}) = 0, the model is not equivalent to a linear regression, as the covariance matrices will still consider species overlap between communities. Therefore, the off-diagnoals in glmm are still not zeroes.
}
\value{
  \item{best_model}{
  INLA output for the best model
  }
  \item{best_model_name}{
  The model identified as the "best".
  }
  \item{predictedfittedresponse}{
  median posterior predictions
  }
  \item{wAIC}{
  wAIC of models (without phylogeny, with the user supplied phylogenetic covariance matrix, and with the optimized phylogeny).
  }
  \item{R2}{
  R-squared calculated as correlations between the observed responses and the mean posterior predictions
  }
  \item{optim_lambda}{
  The value of the optimized \eqn{\lambda}.
  }
  \item{initial_formula}{
  initial formula provided by users.
  }
  \item{without_phylo_model}{
  results on the fixed effects from the model without phylogeny.
  }
  \item{optimized_model}{
  results on the fixed effects from the model after optimizing phylogeny.
  }
  \item{original_VCV_model}{
  results on the fixed effects from the model based on the covariance matrix provided by users.
  }
  \item{Cmatrix}{
  The phylogentic matrix in the optimization model. Note that this is a precision matrix instead of a covariance matrix to fulfill the requirement of INLA
  }
  \item{prediction}{
  predictted values for new data.
  }
}

\author{
Toby P.N. Tsang
}

\seealso{
\code{\link[CPR]{CPR_GLS}}, \code{\link[CPR]{get_comm_pair_r}}, \code{\link[INLA]{inla}}
}

\examples{
library(phytools)
library(tidyverse)
library(CPR)
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

VCV_sp <- vcv(KSR_MLtree)


inla_m <- CPR(formula = bugs~Real.rich+f(comm,model="generic0",Cmatrix=Phylo),
              priors=NULL,
              df = KSR_EF,
              VCV_sp = VCV_sp,
              comm=KSR,
              family="nbinomial",
              optim.lambda = T)

summary(inla_m$best_model)
inla_m$prediction #see predicted relationships at each unique value of the a covariate, while holding others at the mean value.
inla_m$wAIC #wAIC of different models. Note that wAIC_optim has the lowest wAIC.
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
