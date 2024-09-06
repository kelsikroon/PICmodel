
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICmodel

<!-- badges: start -->
<!-- badges: end -->

This `R` package fits prevalence-incidence-cure (PIC) models to
interval-censored data to estimate the cumulative disease risk in a
population with a temporarily elevated risk of disease, e.g., risk of
CIN2+ in HPV-positive women.

In longitudinal screening studies it is possible to observe prevalent,
early or late events during follow-up. In our model, early events are
modelled via a competing risks framework, infections either progress to
the disease state or to a (latent) “cure” state (i.e. viral clearance)
at constant rates. Late events are modelled by adding background risk to
the model. Parameters can depend on individual risk factors and are
estimated with an expectation-maximisation (EM) algorithm with weakly
informative Cauchy priors. More details are given in the accompanying
paper: *eventual link to paper*.

There are five main functions in this package:

- `model.fit`: fits our Prevalence-Incidence-Cure model
- `model.predict`: makes predictions from the model
- `model.simulator`: simulates data under user-specified parameter
  values and covariates
- `score.test.gamma`: performs a score test whether the shape parameter
  of the Gamma distribution for progression is equal to one or not.
- `simulator.gamma`: simulates data where progression follows a Gamma
  rather than exponential distribution and so can have a shape parameter
  not equal to one. Does not allow for covariates, yet.

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
head(sim.dat) # view simulated data
#>        left    right z         age hpv cyt cause      actual
#> 1  0.000000 3.140883 0  0.08416362   0   0     2   0.2436711
#> 2  3.284159 6.302818 0 -0.42056798   1   1     2   3.3496563
#> 3 23.579773      Inf 0  0.73370290   0   0     3 157.8274471
#> 4 24.252462      Inf 0  0.41432009   0   0     3  48.6762930
#> 5 24.170333      Inf 0 -0.68816152   0   1     3 109.1600649
#> 6  0.000000 0.000000 1  0.19925891   0   0     1   0.0000000
sim.fit <- model.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.9698  0.0990 -5.1638 -4.7758
#> g0 intercept   -1.7003  0.0571 -1.8123 -1.5884
#> w0 intercept   -1.2181  0.0682 -1.3518 -1.0844
#> p0 intercept    0.2473  0.0082  0.2312  0.2634
sim.predict <- model.predict(c(), c(), c(), data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)


library(survival) # compare model fit to non-parametric Kaplan-Meier curve 
sim.km.fit <- survfit(Surv(sim.dat$left, sim.dat$right, type='interval2')~1)

# plot PICmodel predictions and KM-curve to compare 
plot(sim.predict[[1]]$Time, sim.predict[[1]]$CR, type='l', xlab='time', ylab='CR', ylim=c(0.2, 0.65))
lines(sim.km.fit$time, 1-sim.km.fit$surv, col='blue')
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />

## Upcoming

- update model predict functions to allow for model fits where the
  general intercept was not estimated.
- update Gamma simulator functions to allow for covariates

## Authors

- **Kelsi Kroon** <k.kroon@amsterdamumc.nl>
