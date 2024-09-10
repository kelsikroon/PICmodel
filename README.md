
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
sim.thetas <- c(-5, -1.6, -1.2, 0.25)
sim.dat <- PICmodel.simulator(3000, c(), c(), c(), sim.thetas, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat) # view simulated data
#>       left    right z      age     age.std hpv cyt cause     actual
#> 1 24.01999      Inf 0 30.85232 -0.82851442   0   1     3  84.454383
#> 2  0.00000 0.000000 1 40.38176 -0.41218222   0   0     1   0.000000
#> 3 24.39051      Inf 0 50.08066  0.01155381   1   1     3 244.513618
#> 4  0.00000 3.023474 0 66.60333  0.73341401   0   0     3   1.271537
#> 5 24.46140      Inf 0 42.73744 -0.30926460   1   0     3  74.034462
#> 6  0.00000 2.882023 0 32.38987 -0.76134047   0   1     2   2.590831

sim.fit <- PICmodel.fit(c(), c(), c(), sim.dat) # fit model to simulated data
sim.fit$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -5.0747  0.0985 -5.2678 -4.8817
#> g0 intercept   -1.5582  0.0543 -1.6646 -1.4517
#> w0 intercept   -1.1330  0.0611 -1.2527 -1.0132
#> p0 intercept    0.2418  0.0081  0.2259  0.2576

sim.predict <- PICmodel.predict(c(), c(), c(), data=sim.dat[1,], time.points = seq(0, 15, 0.5), fit=sim.fit)

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
sim.dat2 <- PICmodel.simulator(5000, c("hpv"), c(), c("cyt"), sim.thetas.cov, show_prob = 0.9, interval=3, include.h=T)
head(sim.dat2) # view simulated data
#>       left    right z      age     age.std hpv cyt cause      actual
#> 1 20.84013      Inf 0 67.80002  0.77756521   0   0     3 180.1536740
#> 2  0.00000 3.079331 0 50.70354  0.04217698   0   0     2   0.5028079
#> 3  0.00000 3.030328 0 35.87145 -0.59581093   1   0     2   1.6694930
#> 4  0.00000 3.758633 0 48.66340 -0.04557788   1   1     2   2.1876578
#> 5  0.00000 0.000000 1 39.34947 -0.44620729   0   1     1   0.0000000
#> 6 23.66946      Inf 0 31.86447 -0.76816762   0   0     3 122.6678531

sim.fit2 <- PICmodel.fit(c("hpv"), c(), c("cyt"), sim.dat2) # fit model to simulated data
sim.fit2$summary # view model fit summary
#>        param theta.hat std.dev   lower   upper
#> h          h   -4.9134  0.0804 -5.0709 -4.7558
#> g0 intercept   -1.5929  0.0487 -1.6884 -1.4975
#> g1       hpv    1.0558  0.0680  0.9225  1.1891
#> w0 intercept   -1.1973  0.0506 -1.2964 -1.0981
#> p0 intercept   -2.9454  0.0887 -3.1193 -2.7716
#> p1       cyt    4.0054  0.1030  3.8035  4.2073

sim.predict2 <- PICmodel.predict(c("hpv"), c(), c("cyt"), data=data.frame(hpv = c(1, 1, 0, 0), cyt=c(1, 0, 1, 0)), 
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
