\name{KSR}
\alias{KSR}

\title{
Data from the Biodiversity-ecosystem functioning experiment in the Koffler Scientific Reserve (KSR)
}
\description{
Species composition and ecosystem functions in each plot. The phylogeny of the fourteen species in the plots is also provided.
}
\usage{
data(KSR)
data(KSR_MLtree)
data(KSR_EF)
}

\format{
  These are linked datasets. \code{KSR} provides species composition data in 88 plots of 14 plant species, while \code{KS_MLtree} describes the phylogenetic relationship between the species. \code{KSR_EF} are the ten ecosystem functions measured in each plot.

  The species in \code{KSR} are \emph{Andropogon gerardii} (ANGE), \emph{Asclepias tuberosa} (ASTU), \emph{Elymus canadensis} (ELCA), \emph{Elymus trachycaulus} (ELTR),
\emph{Desmodium canadense} (DECA), \emph{Lespedeza capitata} (LECA), \emph{Monarda fistulosa} (MOFI), \emph{Penstemon hirsutus} (PEHI), \emph{Pycnanthemum tenuifolium} (PYTE), \emph{Pycnanthemum virginianum} (PYVI), \emph{Rudbeckia hirta} (RUHI), \emph{Schizachyrium scoparium} (SCSC), \emph{Solidago altissima} (SOAL), \emph{Solidago nemoralis} (SONE)

  The meaning of each column in \code{KSR_EF} is explained below:

  \describe{
    \item{\code{Plot}: }{
    Plot identity. This can be ignored. This column also applies to \code{KSR} even though it is absent in the compositional data.}
    \item{\code{Real.rich}: }{
    Plant species richness in the plot.}
    \item{\code{litter2012}: }{
    Amount of litter measured in 2012.}
    \item{\code{ave.biomass}: }{
    Average biomass across 2012-2014.}
    \item{\code{LAI}: }{
    Leaf area index, a simplified dimension of structural complexity.}
    \item{\code{mean.N.change}: }{
    delta 15N change averaged across surface and deep soil.}
    \item{\code{poll_total}: }{
    Total number of pollinators.}
    \item{\code{flwr_total}: }{
    Total number of flowers}
    \item{\code{Mass.loss.2month}: }{
    Decomposition after 2 months.}
    \item{\code{Damage_effect}: }{
    Damage reduction effect.}
    \item{\code{bugs}: }{
    Total number of arthropods.}
    \item{\code{bug.rich}: }{
    Species richness of arthropods.}
    }
    }

\details{
  Note that the compositional data is binary and describe whether the species is planted in each plot initally. The phylogenetic tree is constructed using Maximum Likelihood. See Cadotte et al. (2013) for details.
}

\references{
  Cadotte, M. W. (2013). Experimental evidence that evolutionarily diverse assemblages result in higher productivity. Proceedings of the National Academy of Sciences, 110(22), 8996-9000.
}

\author{
Toby P. N. Tsang
}
