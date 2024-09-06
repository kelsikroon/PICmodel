
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICmodel

<!-- badges: start -->
<!-- badges: end -->

An `R` package for fitting a prevalence-incidence-cure (PIC) model
interval-censored screening data of a population with a temporarily
elevated risk of disease.

The goal of PICmodel is to …

## Installation

You can install the development version of PICmodel from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kelsikroon/PICmodel")
```

## Examples

This is a basic example which shows you how to solve a common problem:

``` r
library(PICmodel)
sim.thetas <- c(-5, -1.6, -1.2, 0.25)
sim.dat <- model.simulator(1000, c(), c(), c(), sim.thetas, show_prob = 0.85, interval=4, include.h=T)
head(sim.dat)
#>        left    right z        age hpv cyt cause      actual
#> 1 23.995091      Inf 0 -0.5468160   0   0     3 124.6424110
#> 2  0.000000 4.454868 0 -0.2464823   1   1     2   0.4868971
#> 3  0.000000 0.000000 1 -0.8121822   0   1     1   0.0000000
#> 4  3.899314 7.833621 0 -0.5162974   0   1     2   4.1474681
#> 5 23.956266      Inf 0  0.8440549   1   0     3 111.5725378
#> 6 23.968028      Inf 0 -0.1025649   0   0     3 168.3403445
sim.fit <- model.fit(c(), c(), c(), sim.dat)
sim.fit$summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.9736  0.1832 -5.3326 -4.6145
#> g0 intercept   -1.5040  0.1141 -1.7276 -1.2803
#> w0 intercept   -1.1706  0.1385 -1.4420 -0.8991
#> p0 intercept    0.2547  0.0147  0.2259  0.2835
```
