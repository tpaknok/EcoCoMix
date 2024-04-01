\name{likelihood.lambda.INLA}
\alias{likelihood.lambda.INLA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizing phylogenetic signals (\eqn{\lambda}) at community levels for INLA models
}
\description{
  An internal function called by CPR. When used with optim, the function compares wAIC of models under different \eqn{\lambda} and identify the best one using L-BFGS-B.
}
\usage{
likelihood.lambda.INLA(inla_formula,
                       data,
                       family = "gaussian",
                       phyV = NULL,
                       comm = NULL,
                       lambda,
                       priors = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inla_formula}{
  an INLA formula
  }
  \item{data}{
  a data frame
}
  \item{family}{
  error distribution (e.g., "gaussian","poisson"). See INLA package for the distribution available.
}
  \item{phyV}{
  a phylogenetic covariance matrix at species levels.
}
  \item{comm}{
  a species-site matrix.
}
  \item{lambda}{
  initial lambda ranging from 0.2 to 0.8/
}
  \item{priors}{
  priors used in the INLA model
}
  \item{\dots}{
  additional arguments passing to optim
}
}
\details{
The initial value is randomized between 0.2 and 0.8, but the results are capped between 0 and 1. The function does not start at the extremes to avoid L-BFGS-B terminating prematurely.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
