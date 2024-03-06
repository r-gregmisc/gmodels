
# gmodels

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/gmodels)](https://CRAN.R-project.org/package=gmodels)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/r-gregmisc/gmodels/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-gregmisc/gmodels/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Tools for fitting linear models that complement those in base `R`. 

Provided functions include:

* `ci` - Compute Confidence Intervals
* `coefFrame` -	Return model parameters in a data frame
* `CrossTable` - Cross Tabulation with Tests for Factor Independence
* `estimable`	- Compute contrasts and estimable linear functions 
* `fast.prcomp` -	Efficient computation of principal components and singular value decomposition
* `fit.contrast` - Compute and test arbitrary contrasts for regression objects
* `glh.test` - Test a General Linear Hypothesis for a Regression Model
* `make.contrasts` - Construct a User-Specified Contrast Matrix

## Installation

Install the released version of gmodels from [CRAN](https://cran.r-project.org) with:

```r
install.packages('gmodels')
```

Install the development version of gmodels from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("r-gregmisc/gmodels")
```

