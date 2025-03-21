#' @importFrom stats D na.omit rexp rnorm runif
#' @importFrom utils tail
#' @importFrom Rlab rbern
NULL

#' Screening data simulator
#'
#' Simulates screening data with user-specified parameters loosely based on the example of cervical cancer screening data. Useful for validating that the model is able to recover the true parameter values. Currently it is only possible to simulate
#' data with the baseline covariates age (continuous between 30 and 70), HPV genotype (HPV16 positive or negative), and cytology (normal/ abnormal).
#'
#' @param n Number of women in the simulated data set.
#' @param prog_model A string in the format "prog ~ x1 + x2 + ...", where x1, x2 are the names of covariates used in the progression rate parameter (must match column name(s) in the input data). Defaults to "prog ~ 1", which is an intercept only model (i.e., without covariates).
#' @param prev_model A string in the format "prev ~ x1 + x2 + ...", where x1, x2 are the names of covariates used in the parameter (probability of prevalent disease) (must match column name(s) in the input data). Defaults to "prev ~ 1", which is an intercept only model (i.e., without covariates).
#' @param params Numerical vector the parameter values to be used in the data simulation (first value is background risk, then l1, l2, pi), including the values of the covariates.
#' @param show_prob A value representing the probability of a woman showing up for the screening visit. Defaults to 0.9.
#' @param interval A value representing the interval between screening rounds (in years). Defaults to 3.
#' @param include.h Indicator variable for whether to background risk in the simulation procedure. Defaults to TRUE.
#' @param covar_settings A vector of settings used to generate covariate values for age, abnormal cytology and HPV16 genotype. They must be supplied in the order age minimum, age maximum, abnormal cytology probability and HPV16 probability. Age is generated uniformly between `age.min` and `age.max`, with default values of 30 and 70 years respectively. Abnormal cytology is generated from a Bernoulli distribution with probability `cyt.prop`, with default probability of 0.4. HPV16 is also generated from a Bernoulli distribution with probability `hpv16.prob', which defaults to 0.3.
#' @return A data frame containing the left and right interval of CIN2+ detection, the indicator of prevalent disease, age (continuous), HPV genotype (HPV16 or other, 1=HPV16), and cytology result (normal or abnormal, 1=abnormal).
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
#' @examples
#' sim.df <- PICmodel.simulator(1000, prog_model = "prog ~ hpv",  c(-5, -1.6, 1, -1.2, 0.25), show_prob = 0.9, interval=3, include.h=T)
#' head(sim.df)
PICmodel.simulator <- function(n, prog_model = "prog ~ 1", prev_model = "prev ~ 1", params, show_prob = 0.9, interval=3, include.h=T, covar_settings = c(age.min = 30, age.max = 70, cyt.prob = 0.4, hpv16.prob = 0.3)){
  age <- runif(n, covar_settings[1], covar_settings[2])
  age.std <- 0.5*(age - mean(age))/sd(age)

  # Cytology Results: this is an indicator variable so 1 means abnormal cytology and 0 means not abnormal (or unknown for z=NA)
  # if they did not show up for screening at time 0 then their cytology result is 0 because it is unknown
  cytology <- rbinom(n, 1, covar_settings[3])

  # HPV genotype (HPV 16 or other) - this is an indicator variable so 1 means they have HPV16 and 0 means other HPV type
  hpv <- rbinom(n, 1, covar_settings[4])

  # model covariates:
  l2_x <- c() # no covariates for clearance/cure

  if (prog_model == 'prog ~ 1'){
    l1_x <- c() # null model without covariates
  }else{
    l1_x <- trimws(strsplit(substr(prog_model, 8, nchar(prog_model)), "[+]")[[1]]) # split model input at every plus sign and ignore the first 7 characters
  }

  if (prev_model == 'prev ~ 1'){
    pi_x <- c() # null model without covariates
  }else{
    pi_x <- trimws(strsplit(substr(prev_model, 8, nchar(prev_model)), "[+]")[[1]]) # split model input at every plus sign and ignore the first 7 characters
  }

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

  covariate_data <- create.covariate.data(data=data.frame(age=age, age.std=age.std, hpv=hpv, cyt=cytology))
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  if (include.h){
    h <- params[1]
    params <- params[-1]
  }

  l1 <- exp(data1 %*% params[1:(n1)])
  l2 <- exp(data2 %*% params[(n1+1):(n1+n2)])
  p <- exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% params[(n1+n2+1):(n1+n2+n3)]))

  #print(p)
  # disease process
  t1 <- ifelse(rbinom(n, 1, p)==1, 0, Inf) # prevalent CIN2/3
  t2 <- ifelse(rbinom(n, 1, l1/(l1+l2))==1, rexp(n, rate=(l1+l2)), Inf) # due to the HPV infection at baseline
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
  screening_times <- data.frame(x1 = ifelse(rbinom(n, 1, show_prob), 0, NA))  # test at 0 years
  for (i in 1:n.tests){
    # loop through number of tests and generate screening time using normal distribution around that interval times test number
    # if interval = 3, then these are generated from beta distributions with means around 3, 6, 9...
    screening_times <- data.frame(cbind(screening_times, xi = ifelse(rbinom(n, 1, show_prob), (i-1)*interval + interval/2 + rbeta(n, 20, 20)*interval, NA) ))
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

  # create left intervals by finding the last value in the list of screens that is smaller than the actual event time
  left <- vapply(screens, function(x) x[Position(function(z) z <= x[length(x)], x[c(1:(length(x))-1)], right=TRUE)], 1)

  # create a list of right intervals by finding the first value of screen times that is greater than the actual event time
  right <- vapply(screens, function(x) x[Position(function(z) z > x[length(x)], x[c(1:(length(x))-1)])], 1)

  z[is.na(left)] <- NA # if left interval is NA then disease was unknown at baseline because it was not checked
  left[is.na(left)] <- 0 # set unknown left intervals to 0 because CIN2/3 could have been there at baseline

  right[left==0 & t==0 & !is.na(z)] <-  0 # if actual time=0 and left=0 and z!=NA ==> disease checked at baseline, then right=0 also
  z[which(right==0)] <- 1 # right is only zero when disease is prevalent (defined above)

  # if the actual time of CIN2/3 development is after the last screening time, then the set the time to Inf
  last_screening_time <- vapply(screens, function(x) tail(x, 2)[1], 1)
  right[screening_times$actual > last_screening_time] <- Inf

  # if right=NA then set to infinity - happens if all screening rounds were NA (very rare, this is to avoid errors in case it happens)
  right[is.na(right)] <- Inf

  return(data.frame(left, right, z = z, age = age, age.std = age.std, hpv = hpv, cyt=cytology, cause=cause, actual=t))
}
