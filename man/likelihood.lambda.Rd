\name{likelihood.lambda}
\alias{likelihood.lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimizing phylogenetic signals (\eqn{\lambda}) at community levels for GLS
}
\description{
  An internal function called by CPR_GLS. When used with optim, the function compares AIC of models under different \eqn{\lambda} and identify the best one using L-BFGS-B.
}
\usage{
  likelihood.lambda(lambda,
  formula,
  df,
  VCV_sp,
  comm,
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
  initial lambda. In CPR_GLS the default is randomly select a value between 0.2-0.8
}
  \item{formula}{
  the formula of the GLS
}
  \item{df}{
  a data frame
}
  \item{VCV_sp}{
  phylogenetic covariance matrix at species levels
}
  \item{comm}{
  a species-site matrix
}
  \item{...}{
  additional arguments passing to optim.
}

}
\details{
The initial value is randomized between 0.2 and 0.8, but the results are capped between 0 and 1. The function does not start at the extremes to avoid L-BFGS-B terminating prematurely.
}
\value{
returning the (\eqn{\lambda}) leading to lowest AIC
}
\author{
Toby P. N. Tsang
}
\seealso{
\code{\link[CPR]{likelihood.lambda}},  \code{\link[CPR]{get_comm_pair_r}, \code{\link[nlme]{gls}}
}
