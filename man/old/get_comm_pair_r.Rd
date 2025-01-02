\name{get_comm_pair_r}
\alias{get_comm_pair_r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_comm_pair_r(comm, VCV_sp, comm_kronecker = NULL, force.PD = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{
%%     ~~Describe \code{comm} here~~
}
  \item{VCV_sp}{
%%     ~~Describe \code{VCV_sp} here~~
}
  \item{comm_kronecker}{
%%     ~~Describe \code{comm_kronecker} here~~
}
  \item{force.PD}{
%%     ~~Describe \code{force.PD} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (comm, VCV_sp, comm_kronecker = NULL, force.PD = F)
{
    require(Matrix)
    require(matrixcalc)
    if (sum(colnames(comm) \%in\% rownames(VCV_sp)) != ncol(comm)) {
        stop("Inconsistent species name between species covariance matrix and community data matrix")
    }
    comm <- comm[, match(colnames(comm), colnames(VCV_sp))]
    if (is.null(comm_kronecker)) {
        product <- as.data.frame(kronecker(as.matrix(comm), as.matrix(comm)))
    }
    else {
        product <- comm_kronecker
    }
    time1 <- Sys.time()
    row_op <- cbind(expand.grid(1:nrow(comm), 1:nrow(comm)))
    row_op$cov <- as.matrix(product) \%*\% c(VCV_sp)
    covM <- matrix(row_op$cov, nrow(comm), nrow(comm))
    corM <- cov2cor(covM)
    time2 <- Sys.time()
    time2 - time1
    if (force.PD == T & is.positive.definite(round(corM, 5)) ==
        F) {
        corM <- as.matrix(nearPD(corM, corr = T, keepDiag = T,
            maxit = 1e+05)$mat)
    }
    result <- list(corM = corM, covM = covM, comm_kronecker = product)
    return(result)
  }
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
