
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PICmodel

<!-- badges: start -->
<!-- badges: end -->

This `R` package fits prevalence-incidence-clearance (PIC) models to
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

- `PICmodel.fit`: fits our Prevalence-Incidence-Cure model
- `PICmodel.predict`: makes predictions from the model
- `PICmodel.simulator`: simulates data under user-specified parameter
  values and covariates
- `score.test.gamma`: performs a score test whether the shape parameter
  of the Gamma distribution for progression is equal to one or not.
- `simulator.gamma`: simulates data where progression follows a Gamma
  rather than exponential distribution and so can have a shape parameter
  not equal to one. Does not allow for covariates, yet.

## Installation

You can install the most recent version of `PICmodel` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kelsikroon/PICmodel")
```

## Examples

This is a basic example in a setting without covariates which
illustrates how to simulate data, fit the model, make predictions from
the plot and compare to a non-parametric cumulative incidence curve:

``` r
library(PICmodel)
sim.thetas <- c(-5, -1.6, -1.2, -3)
sim.dat <- PICmodel.simulator(1000, c(), c(), c(), sim.thetas, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat) # view simulated data
#>        left     right z      age    age.std hpv cyt cause      actual
#> 1 20.958787       Inf 0 52.84289  0.1031094   1   0     3 983.0845280
#> 2  0.000000  6.239495 0 55.11465  0.1997830   1   1     2   3.6913789
#> 3 24.157009       Inf 0 32.87648 -0.7465480   0   1     3  58.2826670
#> 4 24.139272       Inf 0 40.08313 -0.4398735   0   1     3 452.5952674
#> 5  0.000000  3.420154 0 56.78090  0.2706892   0   1     2   0.4380654
#> 6  8.751095 12.161086 0 60.18916  0.4157253   1   0     3   9.3888763

sim.fit <- PICmodel.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.6724  0.1255 -4.9184 -4.4265
#> g0 intercept   -1.5026  0.0829 -1.6651 -1.3401
#> w0 intercept   -1.1412  0.0981 -1.3335 -0.9489
#> p0 intercept   -3.3329  0.1851 -3.6957 -2.9701

sim.predict <- PICmodel.predict(data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)

library(survival) # compare model fit to non-parametric Kaplan-Meier curve 
sim.km.fit <- survfit(Surv(sim.dat$left, sim.dat$right, type='interval2')~1)

# plot PICmodel predictions and KM-curve to compare 
plot(sim.predict[[1]]$Time, sim.predict[[1]]$CR, type='l', xlab='time', ylab='CR')
lines(sim.km.fit$time, 1-sim.km.fit$surv, col='blue')
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />

To add covariates to the model we specify them separately for the
progression, clearance, and prevalence parameters. For example if we
wanted to add HPV16 as a covariate for progression and abnormal cytology
as a covariate for prevalence then we would d the following:

``` r

sim.thetas.cov <- c(-5, -1.6, 1, -1.2, -3, 4)
sim.dat2 <- PICmodel.simulator(1000, c("hpv"), c(), c("cyt"), sim.thetas.cov, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat2) # view simulated data
#>       left   right z      age     age.std hpv cyt cause     actual
#> 1  0.00000 3.29315 0 32.60925 -0.75071887   1   0     2   0.342287
#> 2  0.00000 0.00000 1 50.28429  0.01579363   0   1     1   0.000000
#> 3 24.20113     Inf 0 36.72282 -0.57232613   0   0     3 362.487660
#> 4  0.00000 0.00000 1 30.15812 -0.85701690   0   1     1   0.000000
#> 5  0.00000 0.00000 1 46.04396 -0.16809635   0   1     1   0.000000
#> 6  0.00000 0.00000 1 46.10376 -0.16550329   0   1     1   0.000000

sim.fit2 <- PICmodel.fit(c("hpv"), c(), c("cyt"), sim.dat2) # fit model to simulated data
sim.fit2$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.9112  0.1862 -5.2762 -4.5463
#> g0 intercept   -1.4311  0.1133 -1.6532 -1.2090
#> g1       hpv    0.9781  0.1545  0.6753  1.2809
#> w0 intercept   -1.1209  0.1235 -1.3629 -0.8788
#> p0 intercept   -3.0578  0.2126 -3.4745 -2.6411
#> p1       cyt    4.0741  0.2404  3.6029  4.5453

sim.predict2 <- PICmodel.predict(data=data.frame(hpv = c(1, 1, 0, 0), cyt=c(1, 0, 1, 0)), 
                              time.points = seq(0, 15, 0.5), fit=sim.fit2)

# compare model fit to non-parametric Kaplan-Meier curve 
sim.km.fit1 <- survfit(Surv(left, right, type='interval2')~1, data = sim.dat2[sim.dat2$hpv==1 & sim.dat2$cyt==1,])
sim.km.fit2 <- survfit(Surv(left, right, type='interval2')~1, data = sim.dat2[sim.dat2$hpv==1 & sim.dat2$cyt==0,])
sim.km.fit3 <- survfit(Surv(left, right, type='interval2')~1, data = sim.dat2[sim.dat2$hpv==0 & sim.dat2$cyt==1,])
sim.km.fit4 <- survfit(Surv(left, right, type='interval2')~1, data = sim.dat2[sim.dat2$hpv==0 & sim.dat2$cyt==0,])
```

``` r
# plot PICmodel predictions and KM-curve to compare 
plot(sim.predict2[[1]]$Time, sim.predict2[[1]]$CR, type='l', xlab='time', ylab='CR', ylim=c(0,1), col='darkblue')
lines(sim.km.fit1$time, 1-sim.km.fit1$surv, col='lightblue', lty=2)

lines(sim.predict2[[2]]$Time, sim.predict2[[2]]$CR, col = 'darkgreen')
lines(sim.km.fit2$time, 1-sim.km.fit2$surv, col='darkolivegreen1', lty=2)

lines(sim.predict2[[3]]$Time, sim.predict2[[3]]$CR,col='deeppink2')
lines(sim.km.fit3$time, 1-sim.km.fit3$surv, col='pink1', lty=2)

lines(sim.predict2[[4]]$Time, sim.predict2[[4]]$CR, col='darkorange3')
lines(sim.km.fit4$time, 1-sim.km.fit4$surv, col='orange1', lty=2)

legend("bottomright",
         legend=c("HPV16-pos & abnormal cyt", "HPV16-pos & normal cyt", "HPV16-neg & abnormal cyt", "HPV16-neg & normal cyt"),
         col=c("darkblue", "darkgreen", "deeppink2", "darkorange3"), lty=1, lwd=2, bty = "n", cex=0.75)
```

<img src="man/figures/README-example3-1.png" width="50%" />

## Authors

- **Kelsi R. Kroon** <k.kroon@amsterdamumc.nl>
