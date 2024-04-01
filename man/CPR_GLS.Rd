\name{CPR_GLS}
\alias{CPR_GLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized least squares regressino for community data.
}
\description{
This function conducts GLS for community-level data. This is not recommended because pglmm can better control type-I error, power, and slope estimates. However, for some analyses (e.g., piecewiseSEMs), INLA is not supported yet, so a compromise is to use GLS.
}
}
\usage{
CPR_GLS(formula,
        df,
        VCV_sp,
        comm,
        optim.lambda = TRUE,
        original.VCV = TRUE,
        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a gls formula.
}
  \item{df}{
a data frame.
}
  \item{VCV_sp}{
a covariance matrix at species levels.
}
  \item{comm}{
a species-site matrix.
}
  \item{optim.lambda}{
TRUE/FALSE. If TRUE the optimial lambda and model will also be returned.
}
  \item{original.VCV}{
TRUE/FALSE. If TRUE the model using the supplied VCV_sp will also be returned.
}
  \item{\dots}{
arguments passing to optim.
}
}
\details{
The function calls gls from nlme to conduct the gls, but will automatically calculate the covariance matrix between communities using get_comm_pair_r.
}
\value{
\item {without_phylo_model}{
The model without considering phylogeny (i.e., equivalent to an OLS regression).
}
\item {optimized_model}{
The model with optimized \eqn{\lambda}.
}
\item {original_VCV_model}{
The model based on the phylogenetic covariance matrix provided by users (i.e., without optimization).
\item {AIC}{
AIC of each model.
}
\item {optim_lambda}{
the optimized \eqn{\lambda} (i.e., the \eqn{\lambda} leading to lowest AIC).
}
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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


GLS_lambda_1 <- CPR_GLS(formula = LAI~Real.rich,
                        df = KSR_EF,
                        VCV_sp = VCV_sp,
                        optim.lambda = TRUE,
                        comm = KSR)

GLS_lambda_1$without_phylo_model #models without considering phylogeny
GLS_lambda_1$optimized_model #models without considering phylogeny
GLS_lambda_1$original_VCV_model #models without considering phylogeny
GLS_lambda_1$AIC #AIC of different models
GLS_lambda_1$optim_lambda #optimized lambda if available

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
