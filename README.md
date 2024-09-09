
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
sim.thetas <- c(-5, -1.6, -1.2, 0.25)
sim.dat <- model.simulator(3000, c(), c(), c(), sim.thetas, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat) # view simulated data
#>      left     right  z      age    age.std hpv cyt cause    actual
#> 1  0.0000  2.936714 NA 41.79124 -0.3602243   0   0     3  2.138136
#> 2 21.0141 24.064706  0 44.72637 -0.2331472   1   0     3 22.359147
#> 3  0.0000  0.000000  1 46.90608 -0.1387764   1   0     1  0.000000
#> 4  0.0000  3.167310 NA 64.06899  0.6042943   0   0     1  0.000000
#> 5  0.0000  0.000000  1 54.68161  0.1978662   0   0     1  0.000000
#> 6  0.0000  2.971963  0 46.33303 -0.1635869   0   0     2  1.005333

sim.fit <- model.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.0154  0.1019 -5.2151 -4.8156
#> g0 intercept   -1.6237  0.0535 -1.7285 -1.5188
#> w0 intercept   -1.2448  0.0630 -1.3683 -1.1213
#> p0 intercept    0.2524  0.0082  0.2363  0.2684

sim.predict <- model.predict(c(), c(), c(), data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)

library(survival) # compare model fit to non-parametric Kaplan-Meier curve 
sim.km.fit <- survfit(Surv(sim.dat$left, sim.dat$right, type='interval2')~1)

# plot PICmodel predictions and KM-curve to compare 
plot(sim.predict[[1]]$Time, sim.predict[[1]]$CR, type='l', xlab='time', ylab='CR', ylim=c(0.2, 0.65))
lines(sim.km.fit$time, 1-sim.km.fit$surv, col='blue')
```

<img src="man/figures/README-example-1.png" width="50%" style="display: block; margin: auto;" />

To add covariates to the model we specify them separately for the
progression, clearance, and prevalence parameters. For example if we
wanted to add HPV16 as a covariate for progression and abnormal cytology
as a covariate for prevalence then we would d the following:

``` r

sim.thetas.cov <- c(-5, -1.6, 1, -1.2, -3, 4)
sim.dat2 <- model.simulator(5000, c("hpv"), c(), c("cyt"), sim.thetas.cov, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat2) # view simulated data
#>        left    right z      age      age.std hpv cyt cause     actual
#> 1  0.000000 3.083559 0 38.17603 -0.513049138   1   0     2   1.289715
#> 2 23.455882      Inf 0 49.79450 -0.006752224   0   0     3  72.228146
#> 3  0.000000 0.000000 1 45.22406 -0.205917774   1   0     1   0.000000
#> 4  0.000000 0.000000 1 53.22681  0.142817490   0   1     1   0.000000
#> 5  3.354585 5.988839 0 42.77547 -0.312619582   1   0     2   3.493650
#> 6 23.947902      Inf 0 60.24844  0.448798204   0   0     3 432.020861

sim.fit2 <- model.fit(c("hpv"), c(), c("cyt"), sim.dat2) # fit model to simulated data
sim.fit2$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.8715  0.0769 -5.0222 -4.7208
#> g0 intercept   -1.6057  0.0508 -1.7053 -1.5061
#> g1       hpv    1.0150  0.0681  0.8816  1.1485
#> w0 intercept   -1.1442  0.0525 -1.2471 -1.0413
#> p0 intercept   -2.9472  0.0884 -3.1205 -2.7739
#> p1       cyt    3.8994  0.1019  3.6997  4.0991

sim.predict2 <- model.predict(c("hpv"), c(), c("cyt"), data=data.frame(hpv = c(1, 1, 0, 0), cyt=c(1, 0, 1, 0)), 
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

## Upcoming

- update model predict functions to allow for model fits where the
  general intercept was not estimated.
- update Gamma simulator functions to allow for covariates

## Authors

- **Kelsi Kroon** <k.kroon@amsterdamumc.nl>
