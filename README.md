
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
sim.dat <- PICmodel.simulator(1000, params = sim.thetas, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat) # view simulated data
#>        left     right z      age     age.std hpv cyt cause      actual
#> 1  0.000000  2.919475 0 54.31789  0.17906643   1   1     2   0.2295137
#> 2 24.013075       Inf 0 48.43769 -0.07028581   0   0     3 524.4353027
#> 3 23.518605       Inf 0 41.35735 -0.37053038   1   1     3 498.4729521
#> 4 23.981820       Inf 0 59.89597  0.41560720   1   1     3 452.1208752
#> 5  9.041048 12.277403 0 57.69026  0.32207329   1   0     3   9.2900745
#> 6  0.000000  2.814737 0 39.78798 -0.43708017   0   1     3   1.3926755

sim.fit <- PICmodel.fit(data=sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.9205  0.1469 -5.2084 -4.6325
#> g0 intercept   -1.5805  0.0856 -1.7483 -1.4128
#> w0 intercept   -1.1653  0.1038 -1.3688 -0.9619
#> p0 intercept   -3.5134  0.1972 -3.8999 -3.1269

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
as a covariate for prevalence then we would do the following:

``` r

sim.thetas.cov <- c(-5, -1.6, 1, -1.2, -3, 4)
sim.dat2 <- PICmodel.simulator(1000, prog_model = "prog ~ hpv", prev_model = "prev ~ cyt", sim.thetas.cov, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat2) # view simulated data
#>       left    right z      age     age.std hpv cyt cause     actual
#> 1 23.74135      Inf 0 60.23877  0.43826010   0   1     3  98.454106
#> 2 24.09672      Inf 0 68.94829  0.81634768   0   1     3  62.455445
#> 3  0.00000 2.817943 0 53.68294  0.15366648   1   0     2   1.761897
#> 4 24.04198      Inf 0 51.63299  0.06467654   1   0     3 128.166204
#> 5 24.11061      Inf 0 36.04632 -0.61195336   0   0     3 269.388102
#> 6  0.00000 0.000000 1 62.60351  0.54091556   0   1     1   0.000000

sim.fit2 <- PICmodel.fit(sim.dat2, prog_model = "prog ~ hpv", prev_model = "prev ~ cyt") # fit model to simulated data
sim.fit2$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.0537  0.2185 -5.4819 -4.6254
#> g0 intercept   -1.6742  0.1105 -1.8908 -1.4577
#> g1       hpv    0.9791  0.1510  0.6831  1.2750
#> w0 intercept   -1.2833  0.1271 -1.5325 -1.0342
#> p0 intercept   -2.5920  0.1702 -2.9256 -2.2584
#> p1       cyt    3.4785  0.2039  3.0789  3.8782

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
