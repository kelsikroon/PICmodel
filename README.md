
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
sim.dat <- model.simulator(3000, c(), c(), c(), sim.thetas, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat)
#>       left    right  z        age hpv cyt cause    actual
#> 1 0.000000 0.000000  1  0.3715747   1   0     1 0.0000000
#> 2 0.000000 0.000000  1  0.5315498   0   0     1 0.0000000
#> 3 0.000000 3.084584 NA -0.5111330   0   0     2 0.2785294
#> 4 2.718008 6.431231  0  0.7860201   1   1     2 3.5292963
#> 5 0.000000 3.011932  0 -0.8007818   0   0     2 0.8103534
#> 6 2.958886 9.222339  0 -0.3186591   1   1     2 4.1959996
sim.fit <- model.fit(c(), c(), c(), sim.dat)
sim.fit$summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.0545  0.1051 -5.2605 -4.8485
#> g0 intercept   -1.5275  0.0526 -1.6306 -1.4244
#> w0 intercept   -1.2276  0.0624 -1.3499 -1.1053
#> p0 intercept    0.2516  0.0082  0.2355  0.2676

library(survival)
sim.km.fit <- survfit(Surv(sim.dat$left, sim.dat$right, type='interval2')~1)
sim.predict <- model.predict(c(), c(), c(), data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)
plot(sim.predict[[1]]$Time, sim.predict[[1]]$CR, type='l', xlab='time', ylab='CR', ylim=c(0, 0.6))
lines(sim.km.fit$time, 1-sim.km.fit$surv, col='blue')
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />
