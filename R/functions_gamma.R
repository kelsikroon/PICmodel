#' Gamma Simulator
#'
#' Function to simulate data for a variant of our model that uses the Gamma distribution for Cause 1 (i.e., progression) in the competing risks framework used for incident disease. Note this code is simplified to not include covariates
#' @param n Number of subjects requried in simulated data set
#' @param params Parameter values for simple model without covariates. Note that parameters go in order (1) background risk, (2) lambda_1, (3), lambda_2, and (4) prevalent).
#' @param k Shape parameter for gamma distribution
#' @param interval Screening interval (in years). Defaults to 3.
#' @param include.h Indicator variable for whether to include background risk or not. Defaults to TRUE.
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
simulator.gamma <- function(n, params,  k, interval=3, include.h=T){
  show_prob <- 0.8
  h <- exp(params[1])
  l1 <- exp(params[2])
  l2 <- exp(params[3])
  p <- params[4]
  z <- rep(0, n) # create the indicator variable Z

  # disease process
  t1 <- ifelse(rbern(n, p)==1, 0, Inf) # prevalent CIN2/3
  x <- rgamma(n, shape=k, rate=l1) # simulate from gamma distribution for cause 1 (i.e., progression)
  y <- rexp(n, rate=l2) # simulate from exponential distribution for cause 2 (i.e., clearance)
  t2 <- ifelse(x < y, x, Inf) # take gamma time if it occurs first, otherwise set progression time to infinity

  if (include.h){
    t3 <- rexp(n, rate=h) # due to the background risk (simulate value for everyone)
    t <- pmin(t1, t2, t3) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2, t3), 1, FUN = which.min) # store the cause (just for interest, its never used in the analysis)
  }else{
    t <- pmin(t1, t2) # keep the minimum value as the actual event time
    cause <- apply(data.frame(t1, t2), 1, FUN = which.min)
  }


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
  # create a list of the screening times with NA values removed for each subject
  #screens <- lapply(seq(1, n), function(x) unlist((screening_times)[x,] ))
  screens <- apply(screening_times, 1, function(x) c(na.omit(unlist(x, use.names=FALSE))))
  # create left intervals by finding the last value in the list of screens that is smaller than the actual event time ß
  left <- vapply(screens, function(x) x[Position(function(z) z <= x[length(x)], x[c(1:(length(x))-1)], right=TRUE)], 1)
  z[is.na(left)] <- NA # if left interval is NA then disease was unknown at baseline because it was not checked
  left[is.na(left)] <- 0 # set unknown left intervals to 0 because CIN2/3 could have been there at baseline

  # create a list of right intervals by finding the first value in the list of screen times that is greater than the actual event time
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

  return(data.frame(left, right, z = z,  cause=cause, actual=t))
}



#' Gamma Score test
#'
#' Function to compute p-value of the score test given model fit, step size (dh), order (accuracy, 1, 2 or 4) and data.
#'
#' @param fit Output from model fit.
#' @param dh Step size used in numerical differentiation (e.g., 1e-5, 1e-10)
#' @param accuracy Order of accuracy in numerical differentiation, options are 1, 2, or 4.
#' @param data Data used to fit the model.
#' @return  The output is a list containing the following elements:
#' \itemize{
#' \item grad: value of the gradient at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
#' \item inv.hess: value of the inverse hessian at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
#' \item hess: value of the hessian at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
#' \item loglik: value of the log-likelihood at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
#' \item p.val: p-value of the score test for whether Cause 1 follows a Gamma rather than Exponential distribution (tests whether shape parameter equals 1 or not).
#' }
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
score.test.gamma <- function(fit, dh, accuracy, data){
  #require("stringr")
  fixed.h <- fit$fixed.h
  n <- length(data$right)

  n1 <- length(fit$model[[1]]) +1 # number of covariates for l1 (progression)
  n2 <- length(fit$model[[2]]) +1 # number of covariates for l2 (clearance)
  n3 <- length(fit$model[[3]]) +1 # number of covariates for pi (prevalence probability)

  logl <- function(theta.hat, p){ # loglikelihood function for within score test
    left <- data$left
    right <- data$right
    z <- data$z

    if (br){ # if background risk included
      if(!is.null(fixed.h)){ # if a fixed value of h is supplied
        h <- fixed.h
      }else { # if h was estimated as one of the parameter
        h <- theta.hat[1]
        theta.hat <- theta.hat[-1]
      }
    }else{ # if background risk is not included in the model
      h <- NA
    }

    l1_coef <- theta.hat[1:n1] # store coefficients for l1 (progression)
    l2_coef <- theta.hat[(n1+1):(n1+n2)] # store coefficients for l2 (clearance)
    pi_coef <- theta.hat[(n1+n2+1):(n1+n2+n3)] # store coefficients for pi (prevalence)
    coefs <- list(l1_coef, l2_coef, pi_coef) # make list of coefficients
    lengths <- c(n1, n2, n3) # store number of coefficients for each main parameter

    for (param in 1:3){# loop through parameters and add to a temporary expression which we will take exp() of and store at the end of the loop
      temp_expression <- 0
      for (i in 1:lengths[param]){
        if (i ==1) { # if i=1 then only add intercept
          temp_expression <- temp_expression + coefs[[param]][i]
        }else{ # if i>1 then add coefficient of that parameter multiplied by the data column responding to same coefficient
          temp_expression <- temp_expression + coefs[[param]][i]*data[[fit$model[[param]][i-1]]]
        }
      }

      if (param==1){ # store expression for l1 (progression)
        l1s <- exp(temp_expression)
      }else if (param==2){ # store expression for l2 (clearance)
        l2s <- exp(temp_expression)
      }else if (param ==3 & n3 !=1){ # store expression for pi (prevalence) --> use logit if there are covariates
        prevs <- exp(temp_expression)/(1+exp(temp_expression))
      }else if (param==3 & n3==1){ # store expression for pi (prevalence) --> use raw value if there are NOT covariates
        prevs <- rep(temp_expression, n)
      }
    }

    if (br){# if background risk then add the H(t) = 1-exp(-exp(h)*t) function to the incident function
      f <- function(t, h, l1s, l2s, p) return(ifelse(t==Inf, 1, (1-exp(-exp(h)*t)) + exp(-exp(h)*t)*pgamma(t, shape=p, rate = l1s + l2s)*l1s^p/(l1s + l2s)^p))
    }else { # otherwise only add the "G(t)" function to the incident function
      f <- function(t, h, l1s, l2s, p) return(ifelse(t==Inf, 1, pgamma(t, shape=p, rate = l1s + l2s)*l1s^p/(l1s + l2s)^p))
    }
    f.left <- f(left, h, l1s, l2s, p)
    f.right <- f(right, h, l1s, l2s, p)

    llk <- rep(0, length(right)) # create empty vector to store log-likelihood values
    llk[which(z==1)] <- I(log(prevs))[which(z==1)] # contribution from prevalent subjects
    llk[which(z==0 & right<Inf)] <- I(log((1-prevs)*(f.right - f.left)))[which(z==0 & right <Inf)]
    llk[which(z==0 & right==Inf)] <- I(log((1-prevs)*(1- f.left)))[which(z==0 & right==Inf)]
    llk[which(is.na(z) & right<Inf)] <- I(log(prevs + (1-prevs)*f.right))[which(is.na(z) & right<Inf)] # unknown z so left=0 so f.left=0
    llk[which(is.na(z) & right==Inf)] <- 0 # if z=NA then left=0, so contribution is 0 whe right=Inf also
    return(sum(llk))
  }

  mles <- fit$theta.hat # store maximum likelihood estimates as "mles"
  hess <- fit$hess # store maximum likelihood estimates as "mles"
  hess.names <- c("k", rownames(fit$summary)) # store the row/column names of the new hessian
  br <- ifelse(any(str_detect(hess.names, "h")) | !is.null(fixed.h), T, F) # check if background risk was one of the estimated parameters

  num.mle <- length(mles)+1 # number of MLE is plus 1 because now we add shape parameter
  gamma_gradient <- rep(0, num.mle) # only shape is non-zero because the rest are at the MLE
  gamma_hess <- matrix(0,  nrow=num.mle, ncol=num.mle) # create empty hessian with the right number of rows and columns
  gamma_hess[2:num.mle, 2:num.mle] <- hess # fill "inside hessian" with model output


  if (accuracy==4){
    gamma_gradient[1] <- (-logl(mles, p=1+2*dh) + 8*logl(mles, p=1+dh) - 8*logl(mles, p=1-dh) + logl(mles, p=1-2*dh))/(12*dh) # equation (4) from latex doc

    for (i in 1:num.mle){
      if (i ==1){
        gamma_hess[1,1] <- (-logl(mles, p=1+2*dh) + 16*logl(mles, p=1+dh) - 30*logl(mles, p=1) + 16*logl(mles, p=1-dh) - logl(mles, p=1-2*dh))/(12*dh^2)
      }else{
        dhs <- rep(0, length(mles)) # create empty vector of differences
        dhs[i-1] <- dh # only for the parameter under consideration add the small h
        f_pp <- gamma_hess[1,1] # copied from H[1,1] entry of the Hessian
        # second deriv wrt theta_i (equation 5 but for theta):
        f_xx <- (-logl(mles+2*dhs, p=1) + 16*logl(mles+dhs, p=1) - 30*logl(mles, p=1) + 16*logl(mles-dhs, p=1) - logl(mles-2*dhs, p=1))/(12*dh^2)

        hess_entry <- (16*(logl(mles+dhs, p=1+dh) + logl(mles-dhs, p=1-dh)) - (logl(mles+2*dhs, p=1+2*dh) + logl(mles-2*dhs, p=1-2*dh)) -
                         30*logl(mles, p=1))/(24*dh^2) - f_xx/2 - f_pp/2         # equation (13)
        gamma_hess[1, i] <- gamma_hess[i, 1] <- hess_entry
      }
    }
  }else if (accuracy ==2){
    gamma_gradient[1] <- (-logl(mles, p=1-dh) + logl(mles, p=1+dh))/(2*dh)
    for (i in 1:num.mle){
      if (i ==1){
        gamma_hess[1,1] <- (logl(mles, p=1+dh) - 2*logl(mles, p=1) + logl(mles, p=1-dh))/(dh^2)
      }else{
        dhs <- rep(0, length(mles))
        dhs[i-1] <- dh # only for the parameter under consideration add the small h
        f_pp <- gamma_hess[1,1] # copied from H[1,1] entry of the Hessian
        f_xx <- (logl(mles + dhs, p=1) - 2*logl(mles, p=1) + logl(mles - dhs, p=1))/(dh^2)

        hess_entry <- (logl(mles + dhs, p=1 + dh) - 2*logl(mles, p=1) + logl(mles - dhs, p=1 - dh))/(2*dh^2) - f_xx/2 - f_pp/2
        gamma_hess[1, i] <- gamma_hess[i, 1] <- hess_entry
      }
    }
  }else if (accuracy == 1){
    gamma_gradient[1] <- (-logl(mles, p=1) + logl(mles, p=1+dh))/dh
    for (i in 1:num.mle){
      if (i ==1){
        gamma_hess[1,1] <- (logl(mles, p=1) - 2*logl(mles, p=1+dh)  + logl(mles, p=1 + 2*dh))/(dh^2)
      }else{
        dhs <- rep(0, length(mles))
        dhs[i-1] <- dh # only for the parameter under consideration add the small h
        f_pp <- gamma_hess[1,1] # copied from H[1,1] entry of the Hessian
        f_xx <- (logl(mles, p=1) - 2*logl(mles + dhs, p=1) + logl(mles + 2*dhs, p=1))/(dh^2)
        hess_entry <- (logl(mles + 2*dhs, p=1 + 2*dh) - 2*logl(mles + dhs, p=1 + dh) + logl(mles, p=1))/(2*dh^2) - f_xx/2 - f_pp/2
        gamma_hess[1, i] <- gamma_hess[i, 1] <- hess_entry
      }
    }
  }
  test.stat <-  -1*gamma_gradient[1]^2*solve(gamma_hess)[1] ##1/n * gamma_gradient %*% solve(-gamma_hess/n) %*% gamma_gradient
  p.val <- pchisq(test.stat, df = 1, lower.tail = FALSE)[1]
  colnames(gamma_hess) <- rownames(gamma_hess) <- hess.names
  return(list(grad = gamma_gradient[1], inv.hess = solve(gamma_hess), hess=gamma_hess, loglik = logl(mles, p=1), p.val = p.val))
}
