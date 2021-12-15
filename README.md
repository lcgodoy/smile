<!-- README.md is generated from README.Rmd. Please edit that file -->



# smile

## **S**patial **M**isalignment: **I**nterpolation, **L**inkage, and **E**stimation

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/smile)](https://CRAN.R-project.org/package=smile)
[![R-CMD-check](https://github.com/lcgodoy/spmismo/workflows/R-CMD-check/badge.svg)](https://github.com/lcgodoy/smile/actions)
<!-- badges: end -->

The goal of `smile` is to provide an easy to use `R` package that functions to
estimate, predict and interpolate areal data. For estimation and prediction we
use a model-based approach assuming areal data to be an average of an underlying
continuos Gaussian random field. For interpolation, we use the simple areal
interpolation, which is a nonparametric approach. Also, we provide functions to
quantify the uncertainty associated with the simple areal interpolation when the
variance of the observed variables is known.

The package accompanies a web page (powered by
[pkgdown](https://pkgdown.r-lib.org/)) and 5 vignettes.

### Vignettes

1. [Converting `sf` to `spm`
   objects](https://lcgodoy.me/smile/articles/sf-to-spm.html);
2. [Fitting models and making
   predictions](https://lcgodoy.me/smile/articles/fit-and-pred.html);
3. [Simple Areal Interpolation](https://lcgodoy.me/smile/articles/sai.html);
4. [Theory](https://lcgodoy.me/smile/articles/theory.html);
5. [Spatial covariance functions](https://lcgodoy.me/smile/articles/sp-cov-funs.html);

### Installation

The package has not been submitted to [CRAN](https://CRAN.R-project.org)
yet. Therefore, only the development version from [GitHub](https://github.com/)
is available. If you are interested, you can installed using:
```r
install.packages("smile")
```
