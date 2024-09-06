
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICmodel

<!-- badges: start -->
<!-- badges: end -->

`PICmodel`: A prevalence-incidence-cure (PIC) model for
interval-censored data to estimate the cumulative disease risk in a
population with a temporarily elevated risk of disease, e.g., risk of
CIN2+ in HPV-positive women. In longitudinal screening studies it is
possible to observe prevalent, early or late events during follow-up. In
our model, early events are modelled via a competing risks framework,
infections either progress to the disease state or to a (latent) “cure”
state (i.e. viral clearance) at constant rates. Late events are modelled
by adding background risk to the model. Parameters can depend on
individual risk factors and are estimated with an
expectation-maximisation (EM) algorithm with weakly informative Cauchy
priors. More details are given in the accompanying paper: *eventual link
to paper*.

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
#>       left    right z       age hpv cyt cause   actual
#> 1 0.000000 0.000000 1 0.3316258   0   0     1 0.000000
#> 2 0.000000 3.053300 0 0.1106357   0   1     2 1.133672
#> 3 0.000000 0.000000 1 0.3027752   1   0     1 0.000000
#> 4 2.980117 6.196289 0 0.7201633   1   0     2 2.981977
#> 5 0.000000 2.812348 0 0.6568460   0   1     2 1.903825
#> 6 2.725478 6.163698 0 0.2449390   0   1     2 4.741056
sim.fit <- model.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.0183  0.1080 -5.2299 -4.8066
#> g0 intercept   -1.6763  0.0548 -1.7837 -1.5689
#> w0 intercept   -1.2972  0.0680 -1.4305 -1.1639
#> p0 intercept    0.2478  0.0081  0.2319  0.2637
sim.predict <- model.predict(c(), c(), c(), data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)


library(survival) # compare model fit to non-parametric Kaplan-Meier curve 
sim.km.fit <- survfit(Surv(sim.dat$left, sim.dat$right, type='interval2')~1)

# plot PICmodel predictions and KM-curve to compare 
plot(sim.predict[[1]]$Time, sim.predict[[1]]$CR, type='l', xlab='time', ylab='CR', ylim=c(0.2, 0.65))
lines(sim.km.fit$time, 1-sim.km.fit$surv, col='blue')
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />
