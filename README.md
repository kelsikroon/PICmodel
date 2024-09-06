
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
#>        left    right  z        age hpv cyt cause      actual
#> 1  0.000000 2.916200  0  0.5167249   1   0     2   0.6699919
#> 2 24.322822      Inf  0 -0.4898768   0   1     3 147.5717649
#> 3  0.000000 2.995239 NA  0.6841922   0   1     1   0.0000000
#> 4  0.000000 3.063791  0  0.5575474   1   1     2   1.7127338
#> 5  0.000000 3.109975  0 -0.6926260   0   0     2   0.3726588
#> 6  2.937633 5.706490  0 -0.7476236   0   1     2   4.3717067
sim.fit <- model.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.1045  0.1095 -5.3192 -4.8899
#> g0 intercept   -1.6695  0.0537 -1.7748 -1.5643
#> w0 intercept   -1.2949  0.0642 -1.4207 -1.1691
#> p0 intercept    0.2469  0.0082  0.2309  0.2630
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
- update Gamma simulator functions

## Authors

- **Kelsi Kroon** <k.kroon@amsterdamumc.nl>
