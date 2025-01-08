
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CPR <img src="man/figures/logo.png" align="right" height="150" />

<!-- badges: start -->
<!-- badges: end -->

EcoCoMix addresses phylogenetic non-independence between species in
community analyses.

## Installation

You can install the development version of CPR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tpaknok/EcoCoMix")
```

## Main function

The main function of this package is EcoCoMix, which address
compoisitional autocorrelation in regression. Briefly, the function
automatically constructs the compositional (after considering
phylogenetic covariance between species if provided) matrix between
communities based on their composition and the evolutionary history of
the species.

## Background

In many ecological analyses, we obtained different metrics at community
levels and analyzed them using LM/GLM. Nevertheless, this assumes
statistical independence between all communities. Species however is not
independent, given their shared evolutionary history. Indeed, in many
comparative analyses (e.g., correlating different traits), the shared
evolutionary history between species has to be controlled using
Generalized Least Square Regression or Generalized Linear Mixed Model.

The functions in CPR help users to control phylogenetic non-independence
between communities without intensive coding. Users only need 1) species
composition in each community and 2) a phylogenetic tree including all
species in the community data.
