#' @importFrom stats D na.omit rexp rnorm runif
#' @importFrom utils tail
#' @importFrom Rlab rbern
NULL

#' Screening data simulator
#'
#' Simulates screening data with user-specified parameters based on the example of cervical cancer screening data. Useful for validating that the model is able to recover the true parameter values. Currently it is only possible to simulate
#' data with the baseline covariates age (continuous between 3- and 70), HPV genotype (HPV16 positive or negative), and cytology (normal/ abnormal).
#'
#' @param n Number of women in the simulated data set.
#' @param l1_x A vector containing the names of covariates used in the \ifelse{html}{\out{\eqn{\lambda}<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter. Options are "age", "HPV16" and "cytology".
#' @param l2_x A vector containing the names of covariates used in the \eqn{\lambda_2} (clearance rate) parameter. Options are "age", "HPV16" and "cytology".
#' @param pi_x A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease). Options are "age", "HPV16" and "cytology".
#' @param params Numerical vector the parameter values to be used in the data simulation (first value is background risk, then l1, l2, pi)
#' @param show_prob A value representing the probability of a woman showing up for the screening visit. Defaults to 0.9.
#' @param interval A value representing the interval between screening rounds (in years). Defaults to 3.
#' @param include.h Indicator variable for whether to background risk in the simulation procedure. Defaults to TRUE.
#' @return A data frame containing the left and right interval of CIN2+ detection, the indicator of prevalent disease, age (continuous), HPV genotype (HPV16 or other, 1=HPV16), and cytology result (normal or abnormal, 1=abnormal).
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
#loosely based on cervical cancer data
screening.simulator <- function(n, l1_x, l2_x, pi_x, params, show_prob = 0.9, interval=3, include.h=T){
  age <- runif(n, 30, 70)
  age <- 0.5*(age - mean(age))/sd(age)

  # Cytology Results: this is an indicator variable so 1 means abnormal cytology and 0 means not abnormal (or unknown for z=NA)
  # if they did not show up for screening at time 0 then their cytology result is 0 because it is unknown
  cytology <- rbern(n, 0.4)

  # HPV genotype (HPV 16 or other) - this is an indicator variable so 1 means they have HPV16 and 0 means other HPV type
  hpv <- rbern(n, 0.3)

  create.covariate.data <- function(data){
    data[['intercept']] <- rep(1, length(dim(data)[1]))
    data1 <- as.matrix(data[,c("intercept", l1_x)])
    data2 <- as.matrix(data[,c("intercept", l2_x)])
    data3 <- as.matrix(data[,c("intercept", pi_x)])
    return(list(data1, data2, data3))
  }
  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  covariate_data <- create.covariate.data(data=data.frame(age=age, hpv=hpv, cyt=cytology))
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  if (include.h){
    h <- params[1]
    params <- params[-1]
  }

  l1 <- exp(data1 %*% params[1:(n1)])
  l2 <- exp(data2 %*% params[(n1+1):(n1+n2)])
  if (n3 == 1){
    # if pi has no covariates then the value of p is the same for all women
    p <- params[n1+n2+n3]
  }else if (n3 > 1){
    p <- exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)]))
  }

  # disease process
  t1 <- ifelse(rbern(n, p)==1, 0, Inf) # prevalent CIN2/3
  t2 <- ifelse(rbern(n, l1/(l1+l2))==1, rexp(n, rate=(l1+l2)), Inf) # due to the HPV infection at baseline
  if (include.h) {
    t3 <- rexp(n, rate=exp(h)) # due to the background risk (simulate value for everyone)
    t <- pmin(t1, t2, t3) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2, t3), 1, FUN = which.min)
  }else{
    t <- pmin(t1, t2) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2), 1, FUN = which.min)
  }

  # create observation process
  n.tests <- round(25/interval)  # have at least 25 years of screening
  screening_times <- data.frame(x1 = ifelse(rbern(n, show_prob), 0, NA))  # test at 0 years
  for (i in 1:n.tests){
    # loop through number of tests and generate screening time using normal distribution around that interval times test number
    # if interval = 3, then these are generated from beta distributions with means around 3, 6, 9...
    screening_times <- data.frame(cbind(screening_times, xi = ifelse(rbern(n, show_prob), (i-1)*interval + interval/2 + rbeta(n, 20, 20)*interval, NA) ))
  } #
  #print(head(screening_times))
  screening_times[!is.na(screening_times) & screening_times<0] <- 0


  screening_times$actual <- t
  screening_times$cause <- cause
  screening_times$cause <- NULL
  # create a list of the screening times with NA values removed for each subject
  # the code makes the following change: 0   NA    t2    t3    NA  ---> 0    t2    t3
  # this allows us to find the interval where CIN2/3 was detected
  if (show_prob ==1) {
    screens <- lapply(seq(1, n), function(x) unlist((screening_times)[x,] ))
  }
  else {
    screens <- apply(screening_times, 1, function(x) c(na.omit(unlist(x, use.names=FALSE))))
  }


  z <- rep(0, n) # create the indicator variable Z

  # create left intervals by finding the last value in the list of screens that is smaller than the actual event time ß
  left <- vapply(screens, function(x) x[Position(function(z) z <= x[length(x)], x[c(1:(length(x))-1)], right=TRUE)], 1)

  # if left interval is NA then disease was unknown at baseline because it was not checked
  z[is.na(left)] <- NA

  left[is.na(left)] <- 0 # set unknown left intervals to 0 because CIN2/3 could have been there at baseline

  # create a list of right intervals by finding the first value in the
  # list of screen times that is greater than the actual event time
  right <- vapply(screens, function(x) x[Position(function(z) z > x[length(x)], x[c(1:(length(x))-1)])], 1)

  # if the actual event time t=0 and left interval l=0 and the indicator is not unknown
  # (meaning disease was checked at baseline), then the right interval is also zero
  right[left==0 & t==0 & !is.na(z)] <-  0

  z[which(right==0)] <- 1 # right is only zero when disease is prevalent (defined above)

  # if the actual time of CIN2/3 development is after the last screening time, then the set the time to Inf
  last_screening_time <- vapply(screens, function(x) tail(x, 2)[1], 1)
  right[screening_times$actual > last_screening_time] <- Inf

  # if the right interval is NA then set it to infinity - this happens if all screening
  # rounds were NA (very rare, this is just to avoid errors in case it happens)
  right[is.na(right)] <- Inf

  return(data.frame(left, right, z = z, age = age, hpv = hpv, cyt=cytology, cause=cause, actual=t))
}
