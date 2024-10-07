
#' Fit the prevalence-incidence-cure model
#'
#' This function fits our Prevalence-Incidence-Cure (PIC) mixture model to interval-censored screening data and obtains parameter estimates.
#' It is possible for the user to select the covariates that will be used for each parameter.
#'
#' @param l1_x A vector containing the names of covariates used in the progression rate parameter (must match column name(s) in the input data). Can be left blank for empty model.
#' @param l2_x A vector containing the names of covariates used in the clearance rate parameter (must match column name(s) in the input data). Can be left blank for empty model.
#' @param pi_x A vector containing the names of covariates used in the parameter (probability of prevalent disease) (must match column name(s) in the input data). Can be left blank for empty model.
#' @param data Data used to fit the model containing columns for each term in l1_x, l2_x and pi_x. The first three columns must be: \itemize{
#' \item left interval giving the last time a patient was seen without the disease
#' \item right interval giving the first time a patient was seen with the disease (can be Inf for right-censored cases)
#' \item z indicator for prevalent/incident disease. The data must contain observed prevalent (z=1) and incident (z=0) cases, however some cases may be set to unknown (i.e., z=`NA`) for instance if there were no visits between baseline and disease detection.
#' }
#' @param short.runs Number of runs of the 'Short' EM algorithm used to determine initial values for the EM algorithm. Defaults to 20
#' @param short.iter Number of max iterations used in one run of the 'Short' EM algorithm used to determine initial values for the EM algorithm. Defaults to 10
#' @param short.epsilon Convergence criteria used in the 'Short' EM algorithm to determine initial values. Defaults to 0.1.
#' @param epsilon Convergence criteria for the change in log-likelihood value used for stopping the EM algorithm. Defaults to 1e-08.
#' @param silent Indicator variable for whether to print out each iteration or not. Default value set to FALSE.
#' @param init Option to supply initial values instead of starting from random values in the model fit procedure. Defaults to NULL.
#' @param include.h Indicator variable for whether to include background risk in the model or not. Defaults to TRUE.
#' @param h.method Option for background risk initial values, either 2-step procedure ('2step') or bisection method ('bisection'). Defaults to '2step'.
#' @param two.step.h Indicator variable for whether to perform a two-step procedure to obtain initial values when background risk is included. Defaults to TRUE.
#' @param fixed.h Option to supply a fixed value of background risk, in which case it is treated as fixed and not included in parameter estimation and summary statistics. Needs to be on the log scale, e.g., background risk of 0.01 will use $log(0.01) = 4.6$ as input. Defaults to NULL.
#' @param include.priors Indicator variable for whether to include priors in the estimation procedure. Defaults to TRUE.
#' @param prior.type Choice of prior, either "Cauchy" or "t4" (for a t-distribution prior with 4 degrees of freedom) prior. Defaults to the weakly informative Cauchy prior.
#' @param intercept.prog Indicator variable for whether a general intercept is estimated for the progression parameter. Used when adding country/study specific dummy variables instead. Defaults to TRUE.
#' @param intercept.clear Indicator variable for whether a general intercept is estimated for the clearance parameter. Used when adding country/study specific dummy variables instead. Defaults to TRUE.
#' @param intercept.prev Indicator variable for whether a general intercept is estimated for the prevalence parameter. Used when adding country/study specific dummy variables instead. Defaults to TRUE.
#' @return The output is a list containing the following elements:
#' \itemize{
#' \item model: list containing names of covariates used for each parameter
#' \item initial.values: either the initial values determined by the Short EM algorithm process or the supplied initial values
#' \item fixed.h: the fixed value of background risk if supplied, otherwise NULL
#' \item theta.hat: optimum parameter values estimated by the EM algorithm
#' \item num.iterations: number of iterations until the EM algorithm converged
#' \item log.likelihood: value of the log.likelihood at the maximum likelihood estimates
#' \item hess: hessian matrix at the maximum likelihood estimates
#' \item grad: value of the gradient vector at the the maximum likelihood estimates
#' \item std.dev: standard deviation of parameter estimates
#' \item summary: data frame with estimate, std.dev, and 95\% CI for each parameter (useful for data set comparisons)
#' }
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
#'
#' @examples
#' fit <- PICmodel.fit(c(), c(), c(), sim.dat)
#' fit$summary
PICmodel.fit <- function(l1_x = c(), l2_x= c(), pi_x=c(), data, epsilon=1e-08, short.epsilon=1e-1, short.iter=10, short.runs=20, silent=T,  init=NULL,
                         include.h=T, h.method='2step', two.step.h=T, include.priors=T, prior.type = 'cauchy', fixed.h=NULL, intercept.prog = T, intercept.clear = T, intercept.prev=T){
  if (any(!c(l1_x, l2_x, pi_x) %in% colnames(data))) stop("Covariate is not a column name in input data set. Please check")
  sapply(c("Rlab", "dplyr"), require, character.only = TRUE)
  # g(): the competing risks function used in the 'incidence' part of the mixture model
  # Input:
  #   - l1: the rate parameter for progression to CIN2+
  #   - l2: the rate parameter for viral clearance
  #   - t: time
  # Output:
  #   - the cumulative risk of CIN2+ among HPV-positive women wihtout CIN2+ at baseline
  g <- function(l1, l2, t){
    return(l1 / (l1 + l2) * (1 - exp(-(l1 + l2)* t)))
  }

  # create.par(): function to create the parameter names which will be used in the mstep expressions
  # Input:
  #   - n1: number of covariates used for the progression rate parameter + 1
  #   - n2: number of covariates used for the clearance rate parameter + 1
  #   - n3: number of covariates used in the parameter for the probability of prevalent disease + 1
  # Output:
  #   - vector of names of the parameters with (1) g0, g1... for l1, (2) w0, w1... for l2 and (3) p0, p1... for pi (depending on whether there are intercepts or not)
  create.par <- function(n1, n2, n3) {
    if (intercept.prog) par1 <- paste0(c("g"), seq(0, n1 - 1)) else par1 <- paste0(c("g"), seq(1, n1))
    if (intercept.clear) par2 <- paste0(c("w"), seq(0, n2 - 1)) else par2 <- paste0(c("w"), seq(1, n2))
    if (intercept.prev) par3 <- paste0(c("p"), seq(0, n3 - 1)) else par3 <- paste0(c("p"), seq(1, n3))
    return(c(par1, par2, par3))
  }

  # create.covariate.data(): creates matrices for each parameter of observed independent variables (like the 'X' matrices in a simple regression model)
  # Input:
  #   - data: first three columns must be (1) left interval, (2) right interval and (3) z indicator for prevalent/incident disease,
  #           following columns must match the column names given in l1_x, l2_x and pi_x
  # Output:
  #   - list of 3 matrices corresponding to relevant covariates for each of the parameters lambda_1, lambda_2 and pi
  create.covariate.data <- function(data){
    data[['intercept']] <- rep(1, length(dim(data)[1]))
    params.prog <- l1_x
    params.clear <- l2_x
    params.prev <- pi_x

    if(intercept.prog) params.prog <- c("intercept", params.prog)
    if(intercept.clear) params.clear <- c("intercept", params.clear)
    if(intercept.prev) params.prev <- c("intercept", params.prev)

    data1 <- as.matrix(data[,params.prog])
    data2 <- as.matrix(data[,params.clear])
    data3 <- as.matrix(data[,params.prev])

    return(list(data1, data2, data3))
  }

  # log.likelihood(): function for the log-likelihood
  # Input:
  #   - current_par: the current value of the parameters
  #   - data: data frame with columns (1) left, (2) right, and (3) z, representing the left and right interval
  #               and the prevalent/incident indicator variable
  # Output:
  #   - value of the log-likelihood for the current parameter values
  log.likelihood.h <- function(current_par, data, include.h, est.h=T, fixed.h=fixed.h){
    z <- data[[3]] # indicator variable for prevalent or incident disease
    if (include.h){
      if (est.h) {
        h <- exp(current_par[1])
        current_par <- current_par[-1]
      }else {
        h <- exp(fixed.h)
      }
    }else {
      h <-0
    }

    # multiply data by corresponding parameters to calculate l1, l2 and pi
    l1 <- exp((data1) %*% current_par[1:(n1)]) # create current lambda 1 vector for each value in the data
    l2 <- exp((data2) %*% current_par[(n1+1):(n1+n2)]) # create current lambda 1 vector for each value in the data
    if (n3 == 1 & intercept.prev){
      # if there are no covariates for pi (n==1 since only intercept), then calculate pi for each data value
      p <- (data3) %*% current_par[n1+n2+n3]
    }else {
      # if there are covariates for pi then calculate pi for each data value using transformation from logit
      p <- exp((data3) %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp((data3) %*% current_par[(n1+n2+1):(n1+n2+n3)]))
    }

    f <- function(t) return(ifelse(t==Inf, 1, (1-exp(-h*t)) + exp(-h*t)*(l1/(l1 + l2))*(1-exp(-(l1 + l2)*t))))

    llk <- rep(0, length(right)) # create empty vector to store log-likelihood values
    llk[which(z==1)] <- I(log(p))[which(z==1)]
    llk[which(z==0 & right<Inf)] <- I(log((1-p)*(f(right) - f(left))))[which(z==0 & right <Inf)]
    llk[which(z==0 & right==Inf)] <- I(log((1-p)*(1- f(left))))[which(z==0 & right==Inf)]
    llk[which(is.na(z) & right<Inf)] <- I(log(p + (1-p)*f(right)))[which(is.na(z) & right<Inf)]
    llk[which(is.na(z) & right==Inf)] <- 0
    return(sum(llk))
  }

  estep.h <- function(current_par, data, include.h, est.h=T, fixed.h=fixed.h){
    z <- data[[3]] # indicator variable for prevalent or incident disease
    if (include.h){
      if (est.h) {
        h <- exp(current_par[1])
        current_par <- current_par[-1]
      }else {
        h <- exp(fixed.h)
      }
    }else {
      h <-0
    }

    # multiply data by corresponding parameters to calculate l1, l2 and pi
    l1 <- exp((data1) %*% current_par[1:(n1)])
    l2 <- exp((data2) %*% current_par[(n1+1):(n1+n2)])
    if (n3 == 1 & intercept.prev){
      # if pi has no covariates then the value of p is the same for all women
      p <- (data3) %*% current_par[n1+n2+n3]
    }else {
      p <- exp((data3) %*% current_par[(n1+n2+1):(n1+n2+n3)])/(1+exp((data3) %*% current_par[(n1+n2+1):(n1+n2+n3)]))
    }
    # the expected value of z only depends on the right interval since the left is zero when z is unknown
    est <- p/(p + (1 - p)*(1-exp(-h*right) + exp(-h*right)*g(l1, l2, right)))
    est[which(right==Inf)] <- p[which(right==Inf)]
    # the expected value of the latent variable is the estimate above if z is unknown otherwise it is the true value of z
    z[is.na(z)] <-  est[is.na(z)]
    return(z)
  }

  build.expr.h <- function(model_par, include.h, est.h=T, fixed.h=fixed.h){
    # function to 'build' the expected log-likelihood function as an R 'expression' which
    # allows it to be differentiated w.r.t. different parameters in the mstep() function

    # we loop through the names of covariate parameters (e.g., g0, g1, w0, p0, p1, p2) for each of l1, l2 and pi
    # and paste them together with the corresponding data matrix, this string is then used to create the
    # expression which will be differentiated w.r.t each parameter in the mstep() function
    if(include.h & est.h) model_par <- model_par[-1] # remove h as first parameter if h is included as one of the parameters estimated

    # expression for lambda_1 (progression rate parameter)
    expr1 <- ifelse(intercept.prog, model_par[1], "")
    if (length(l1_x) > 0){
      if (intercept.prog) expr1 <- paste0(expr1, "+")
      for (i in 1:length(l1_x)){
        expr1 <- paste0(expr1, ifelse(i!=1, "+", ""), model_par[i + ifelse(intercept.prog, 1, 0)], "*", "data1", i)
    }}

    # expression for lambda_2 (clearance rate parameter)
    expr2 <- ifelse(intercept.clear, model_par[n1+1], "")
    if (length(l2_x) >0){
      if (intercept.clear) expr2 <- paste0(expr2, "+")
      for (i in 1:length(l2_x)){
        expr2 <- paste0(expr2, ifelse(i!=1, "+", ""), model_par[n1+i+ifelse(intercept.clear, 1, 0)], "*", "data2", i)
    }}

    # expression for pi (prevalent probability parameter)
    expr3 <- ifelse(intercept.prev, model_par[n1+n2+1], "")
    if (length(pi_x) > 0){
      if (intercept.prev)  expr3 <- paste0(expr3, "+")
      for (i in 1:length(pi_x)){
        expr3 <- paste0(expr3, ifelse(i!=1, "+", ""), model_par[n1+n2+i+ifelse(intercept.prev, 1, 0)], "*", "data3", i)
    }}

    p.expr <- ifelse(n3==1 & intercept.prev, str2expression(paste0("z*log(", expr3, ") + (1-z)*log(1-", expr3, ")")),str2expression(paste0("z*(", expr3, ") - log(1 + exp(",expr3, "))")))
    f.expr <- str2expression(paste0("(1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",expr1 ,")/(exp(", expr1, ") + exp(",expr2,")))*(1-exp(-(exp(",expr1,") + exp(", expr2,"))*t))"))
    return(list(p.expr, f.expr))
  }


  mstep.h <- function(current_par, expected_z, data, include.h, est.h=T, fixed.h=fixed.h){
    z <- expected_z # expected value of z = output from the estep.h() function
    m1 <- n1+n2
    if (include.h){
      if (est.h) {
        m1 <- n1+n2 + 1 # if estimating h then add one to number of parameters
      }else {
        h <- fixed.h
      }
    }else {
      h <- -Inf # if h is not estimated set at zero
    }

    # assign the values of current_par to vectors with names from the 'pars' list (this will or won't have h and the number of parameters will correspond)
    for (i in 1:length(pars)){
      assign(pars[i], current_par[i], envir=environment())
    }

    f <- function(t) return(ifelse(t==Inf, 1, (eval(f.expr))))     # incident function

    deriv.f <- function(t, param, order=1) { # derivative of incident function
      if (order ==1){
        return(ifelse(t==Inf, 0, (eval(D(f.expr, param)))))
      }else if (order==2){
        return(ifelse(t==Inf, 0, (eval(D(D(f.expr, param[1]), param[2])))))
      }
    }

    model.exprs <- build.expr.h(pars, include.h, est.h)
    p.expr <- model.exprs[[1]] # prevalent function
    f.expr <- model.exprs[[2]] # incident function

    deriv.prior <- function(param, order=1){
      if(! include.priors) return(0)

      if (prior.type == 'cauchy'){ # cauchy prior on the parameter: param = pars[i], mean=0, t=2.5
        if (order==1){ # derivative of log cauchy prior
          #cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param, "^2 + 2.5^2))"))
          cauchy.prior <- str2expression(paste0("1/(pi*2.5*(1+ (", param ,"/2.5)^2 ))"))
          return(eval(D(cauchy.prior, param)))
        }else if (order==2){
          #cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param[1], "^2 + 2.5^2))"))
          cauchy.prior <- str2expression(paste0("1/(pi*2.5*(1+ (", param[1],"/2.5)^2 ))"))
          return(eval(D(D(cauchy.prior, param[1]), param[2])))
        }
      }else if (prior.type == 't4'){
        if (order==1){ # derivative of log t_4 prior
          # str2expression(paste0("(-6*", param, ")/(5 + ", param, "^2)"))
          t4.prior <- str2expression(paste0("3/(8*(1 + 1/4 * ", param,"^2)^(5/2))"))
          return(eval(D(t4.prior, param)))
        }else if (order==2){
          # str2expression(paste0("(-6*(-", param[1], "^2 + 5)/( 5 + ", param[1], "^2)^2"))
          t4.prior <- str2expression(paste0("3/(8*(1 + 1/4 * ", param[1],"^2)^(5/2))"))
          return(eval(D(D(t4.prior, param[1]), param[2])))
        }

      }
    }

    hess <- matrix(rep(0, (length(current_par))^2), nrow=length(current_par)) # empty hessian matrix
    grad <- rep(0, length(current_par)) # empty gradient matrix

    for (i in 1:m1){
      # we loop through the parameter list and calculate the first derivatives for gradient vector for the l1 and l2 parameters (and h if estimated)
      grad[i] <- sum(((1-z)*(deriv.f(right, pars[i]) - deriv.f(left, pars[i]))/(f(right) - f(left)))[z!=1]) + deriv.prior(pars[i])
      for (j in 1:i){
        # calculate 2nd derivatives for the Hessian matrix for h, l1 and l2 parameters;
        # the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(((1-z)*(f(right) - f(left))^(-2)*
                                         ((deriv.f(right, c(pars[i], pars[j]), 2) - deriv.f(left, c(pars[i], pars[j]),2))*(f(right) - f(left)) -
                                            (deriv.f(right, pars[i], 1) - deriv.f(left, pars[i],1))*(deriv.f(right, pars[j],1) - deriv.f(left, pars[j], 1))))[z!=1]) +
          deriv.prior(c(pars[i], pars[j]), order=2)

      }
    }
    # parameters for pi (for example p0, p1, p2, p3, p4) do not depend on left/right so we can calculate it separately
    for (i in (m1+1):length(current_par)){
      # calculate first derivatives for gradient vector for the pi parameters
      grad[i] <- sum(eval(D(p.expr, pars[i]))) + deriv.prior(pars[i])
      for (j in (m1+1):i){
        # calculate 2nd derivatives for Hessian matrix for the pi parameters; the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(eval(D(D(p.expr, pars[i]), pars[j]))) + deriv.prior(c(pars[i], pars[j]), order=2)
      }
    }
    new_par <-  current_par - solve(hess)%*%grad # single Newton step to calculate updated parameters
    return(list(as.vector(new_par), hess, grad)) # return the new parameters
  }

  # numerical.hess(): function to calculate the hessian using numerical methods
  numerical.hess <- function(mles, data, include.h, est.h=T, fixed.h=fixed.h){
    for (i in 1:length(pars)){    # assign the values of current_par to vectors with names from the 'pars' list (this will or won't have h and the number of parameters will correspond)
      assign(pars[i], mles[i], envir=environment())
    }

    deriv.prior <- function(param, order=1){
      if(! include.priors) return(0)

      if (prior.type == 'cauchy'){ # cauchy prior on the parameter: param = pars[i], mean=0, t=2.5
        if (order==1){ # derivative of log cauchy prior
          #cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param, "^2 + 2.5^2))"))
          cauchy.prior <- str2expression(paste0("1/(pi*2.5*(1+ (", param ,"/2.5)^2 ))"))

          return(eval(D(cauchy.prior, param)))
        }else if (order==2){
          #cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param[1], "^2 + 2.5^2))"))
          cauchy.prior <- str2expression(paste0("1/(pi*2.5*(1+ (", param[1],"/2.5)^2 ))"))
          return(eval(D(D(cauchy.prior, param[1]), param[2])))
        }
      }else if (prior.type == 't4'){
        if (order==1){ # derivative of log t_4 prior
          # str2expression(paste0("(-6*", param, ")/(5 + ", param, "^2)"))
          t4.prior <- str2expression(paste0("3/(8*(1 + 1/4 * ", param,"^2)^(5/2))"))
          return(eval(D(t4.prior, param)))
          #return(eval(str2expression(paste0("(-6*", param, ")/(5 + ", param, "^2)"))))
        }else if (order==2){
          # str2expression(paste0("(-6*(-", param[1], "^2 + 5)/( 5 + ", param[1], "^2)^2"))
          t4.prior <- str2expression(paste0("3/(8*(1 + 1/4 * ", param[1],"^2)^(5/2))"))
          return(eval(D(D(t4.prior, param[1]), param[2])))
          #return(eval(str2expression(paste0("(-6*(-", param[1], "^2 + 5)/( 5 + ", param[1], "^2)^2"))))
        }

      }
    }

    dh <- 1e-4 # step size for numerical differentiation
    hess.size <- length(mles)
    hess <- matrix(NA, nrow=hess.size, ncol= hess.size) # initialise empty hessian matrix
    for (i in 1:hess.size){
      for (j in 1:hess.size){
        dh_i <- dh_j <- dh_ij <- rep(0, hess.size) # create vector of zeros the same length as number of parameters
        dh_i[i] <- dh # only for the parameter i under consideration add the small h
        dh_j[j] <- dh # only for the parameter j under consideration add the small h
        dh_ij[c(i,j)] <- dh # for both i and j add the small h

        if (i ==j){ # only the diagonal entries have non-zero prior derivatives
          hess[i,i] <- (-log.likelihood.h(mles+2*dh_i,  data, include.h, est.h, fixed.h) + 16*log.likelihood.h(mles+dh_i,  data, include.h, est.h, fixed.h) - 30*log.likelihood.h(mles,  data, include.h, est.h, fixed.h) +
                          16*log.likelihood.h(mles-dh_i,  data, include.h, est.h, fixed.h) - log.likelihood.h(mles-2*dh_i,  data, include.h, est.h, fixed.h))/(12*dh^2) + deriv.prior(c(pars[i], pars[i]), order=2)

        }else{ # accuracy to the order 4
          f_ii <- (-log.likelihood.h(mles+2*dh_i,  data, include.h, est.h, fixed.h) + 16*log.likelihood.h(mles+dh_i,  data, include.h, est.h, fixed.h) - 30*log.likelihood.h(mles,  data, include.h, est.h, fixed.h) +
                     16*log.likelihood.h(mles-dh_i,  data, include.h, est.h, fixed.h) - log.likelihood.h(mles-2*dh_i,  data, include.h, est.h, fixed.h))/(12*dh^2) + deriv.prior(c(pars[i], pars[i]), order=2)
          f_jj <- (-log.likelihood.h(mles+2*dh_j,  data, include.h, est.h, fixed.h) + 16*log.likelihood.h(mles+dh_j,  data, include.h, est.h, fixed.h) - 30*log.likelihood.h(mles,  data, include.h, est.h, fixed.h) +
                     16*log.likelihood.h(mles-dh_j,  data, include.h, est.h, fixed.h) - log.likelihood.h(mles-2*dh_j,  data, include.h, est.h, fixed.h))/(12*dh^2) + deriv.prior(c(pars[j], pars[j]), order=2)

          hess_entry <- (16*(log.likelihood.h(mles+dh_ij,  data, include.h, est.h, fixed.h) + log.likelihood.h(mles-dh_ij,  data, include.h, est.h, fixed.h)) -
                           (log.likelihood.h(mles+2*dh_ij,  data, include.h, est.h, fixed.h) + log.likelihood.h(mles-2*dh_ij,  data, include.h, est.h, fixed.h)) -
                           30*log.likelihood.h(mles,  data, include.h, est.h, fixed.h))/(24*dh^2) - f_ii/2 - f_jj/2

          hess[i, j] <- hess[j, i] <- hess_entry
        }
      }
    }
    return(hess)
  }

  # em.function.h(): combining the E- and M-step and repeating until the parameter estimates converge
  em.function.h <- function(init, data, include.h, est.h=T, fixed.h=fixed.h){
    new_theta <- init
    old_llk <- 0
    new_llk <- 100 # make big to start with
    iter <- 1

    while ((abs(new_llk - old_llk) > epsilon)){
      current_theta <- new_theta
      old_llk <- new_llk
      next_em_step <- mstep.h(current_theta, estep.h(current_theta, data, include.h, est.h, fixed.h=fixed.h), data, include.h, est.h, fixed.h=fixed.h)
      new_theta <- as.vector(next_em_step[[1]])
      new_llk <- log.likelihood.h(new_theta, data, include.h, est.h, fixed.h=fixed.h)
      if (!silent & (iter %/% 10 ==0 | iter %% 10 ==0)) print(round(c(new_theta, new_llk), 4))
      iter <- iter + 1
      if (iter > 1000){
        break
      }
    }
    hess <- numerical.hess(new_theta, data, include.h, est.h, fixed.h)
    names(new_theta) <- pars
    std.dev <- round(sqrt(diag(solve(-hess))),4)
    param <- c(if(include.h & est.h){"h"}, if(intercept.prog) {"intercept"}, l1_x, if(intercept.clear) {"intercept"}, l2_x, if(intercept.prev) {"intercept"}, pi_x)
    summary <- round(data.frame(theta.hat = new_theta, std.dev, lower = new_theta - 1.96*std.dev, upper = new_theta + 1.96*std.dev), 4)
    summary <- cbind(param, summary)
    rownames(summary) <- names(new_theta)
    return(list(model = list(l1_x, l2_x, pi_x), initial.values = init, fixed.h = fixed.h, theta.hat = new_theta, num.iterations=iter,
                log.likelihood = new_llk, hess=hess, grad = next_em_step[[3]], std.dev=std.dev, summary=summary))
  }

  # short.em():  lots of short runs of the em function to generate appropriate starting values
  short.em <- function(data, init.h.only, include.h, init, est.h=T, fixed.h){
    if(!init.h.only){ #initial values for only l1, l2, p
      new_theta <- log(runif(n1+n2+n3))
      if (n3==1 & intercept.prev) new_theta[n1+n2+n3] <- exp(new_theta[n1+n2+n3]) # if only one parameter for prevalence, then no covariates ==> they're not on log scale
      if(include.h & est.h) { # only generate if estimating h
        new_theta <- c(log(0.001), new_theta) #c(log(runif(1, 0, 0.02)), new_theta)
      }
    }else if (!is.null(fixed.h)){
      new_theta <- c(fixed.h, init)
    }else{# initial values for background risk using previously found values of l1, l2, p
      new_theta <- c(log(0.001), init) # c(log(runif(1, 0, 0.02)), init)
    }
    new_llk <- 100
    old_llk <- 0
    iter <- 1
    while ((abs(new_llk - old_llk) > short.epsilon) & iter < short.iter){
      current_theta <- new_theta
      old_llk <- new_llk
      iter <- iter + 1
      new_theta <- try(as.vector(mstep.h(current_theta, estep.h(current_theta, data, include.h, est.h, fixed.h), data, include.h, est.h, fixed.h)[[1]]), silent=T)
      if (class(new_theta) == 'try-error') return(short.em(data, init.h.only, include.h, init, est.h, fixed.h))
      new_llk <- log.likelihood.h(new_theta, data, include.h, est.h, fixed.h)
      if (abs(new_llk) == Inf | is.nan(new_llk) | is.na(new_llk)) return(short.em(data, init.h.only, include.h, init, est.h, fixed.h))
    }
    return(list(theta.hat = new_theta, log.likelihood = new_llk))
  }

  # init.generator(): performs short runs of the em function using short.em and then returns the most common set of
  # starting values with the highest loglikelihood as the initial values
  init.generator <- function(data, init.h.only, include.h, init, est.h=T, fixed.h){
    if (!silent) pb <- txtProgressBar(min = 0, max = short.runs, style = 3, width = 50, char = "=")
    short.inits <- list()
    for(k in 1:short.runs) {
      short.inits[[k]] <- short.em(data, init.h.only, include.h, init, est.h, fixed.h)[c("theta.hat", "log.likelihood")]
      if (!silent) setTxtProgressBar(pb, k)
    }
    if (!silent) close(pb)
    short.inits.mat <- matrix(unlist(short.inits), nrow=length(short.inits), byrow=T)
    ncols <- ncol(short.inits.mat)
    # find the set of parameter values that results in the maximum likelihood
    init.counts <- cbind(round(short.inits.mat[rev(order(short.inits.mat[,ncols])),1:(ncols-1)],1),
                         round(short.inits.mat[rev(order(short.inits.mat[,ncols])),ncols])) %>%
      data.frame()%>% group_by_all() %>% dplyr::count() %>% as.data.frame()

    if (all(init.counts$n==1)){
      # if there are no repeated starting values, we take the one with the highest likelihood
      inits <- short.inits.mat[which.max(short.inits.mat[,dim(short.inits.mat)[2]]),1:(dim(short.inits.mat)[2]-1)]
    }else{
      # of the starting values which are repeated more than once, take the one with the highest likelihood
      inits <- init.counts %>% arrange(., n) %>% filter(n > 1) %>% slice_max(n=1, order_by = get(noquote(paste0("X", ncols)))) %>% head(1) %>% as.numeric() %>% head(-2)
    }
    if (any(inits ==0)){
      inits[inits==0] <- 0.001
    }
    return(inits)
  }

  left <- data[[1]]
  right <- data[[2]] # right intervals from the data

  n1 <- length(l1_x) + ifelse(intercept.prog, 1, 0) # number of parameters for lambda_1 (add 1 for intercept if estimated)
  n2 <- length(l2_x) + ifelse(intercept.clear, 1, 0) # number of parameters for lambda_2 (add 1 for intercept if estimated)
  n3 <- length(pi_x) + ifelse(intercept.prev, 1, 0) # number of parameters for pi (add 1 for intercept if estimated)

  covariate_data <- create.covariate.data(data)
  data1 <- covariate_data[[1]]
  data2 <- covariate_data[[2]]
  data3 <- covariate_data[[3]]

  {
  if ((n1 !=1 & intercept.prog) | (n1 > 0 & !intercept.prog)){ # if there are covariates for progression
    for (i in 1:length(l1_x)){
      assign(paste0("data1", i), data1[,i + ifelse(intercept.prog, 1, 0)], envir=environment())
    }}
  if ((n2 !=1 & intercept.clear) | (n2 > 0 & !intercept.clear)){ #if there are covariates for clearance
    for (i in 1:length(l2_x)){
      assign(paste0("data2", i), data2[,i + ifelse(intercept.clear, 1, 0)], envir=environment())
    }}
  if ((n3 !=1 & intercept.prev) | (n3 > 0 & !intercept.prev)){ # if there are covariates for prevalence
    for (i in 1:length(pi_x)){
      assign(paste0("data3", i), data3[,i + ifelse(intercept.prev, 1, 0)], envir=environment())
    }}

  pars <- create.par(n1, n2, n3) # creates a vector of the names of parameters that need to be estimated
  }

  if (!is.null(fixed.h)) {
    if (h.method != '2step') two.step.h <- F
    include.h <- T
    est.h <- F
  }else if(include.h){
    if (h.method == '2step') two.step.h <- T
    else if (h.method == 'bisection') two.step.h <- F
    est.h <- T
  }else{
    est.h <- F
  }

  df <- function(f, x, ...) {
    # calculate the derivative of the log-likelihood using the formal definition of the derivative
    dx <- 0.00001
    d <- (f(fixed.h = x + dx, ...) - f(fixed.h=x-dx, ...)) / 2*dx
    return(d)
  }


  bisection.em <- function(a, b, lower.fit, upper.fit){
    max.iter <- 50
    iter <- 1
    tol <- 1e-5
    mid.init <- lower.fit$theta.hat
    while (abs(a - b) > tol && max.iter > iter){
      c <- (a+b)/2

      mid.fit <- PICmodel.fit(l1_x , l2_x, pi_x, data, epsilon, short.epsilon, short.iter=5, short.runs=5, silent=T,  init = mid.init,
                              include.h, h.method="", two.step.h=F, include.priors, prior.type, fixed.h=log(c), intercept.prog, intercept.clear, intercept.prev)

      lower.res <- df(log.likelihood.h, x = log(a), current_par = lower.fit$theta.hat, data=data, include.h = include.h, est.h =F)
      upper.res <- df(log.likelihood.h, x = log(b), current_par = upper.fit$theta.hat, data=data, include.h = include.h, est.h =F)
      mid.res <- df(log.likelihood.h, x = log(c), current_par = mid.fit$theta.hat, data=data, include.h = include.h, est.h =F)

      if (mid.res*upper.res >= 0 ){ # if middle and upper are same sign then upper bound = middle bound ==> the change is between lower and middle
        b <- c # b is the upper bound value
        upper.fit <- mid.fit
        mid.init <- mid.fit$theta.hat
      }else if (mid.res*lower.res >= 0){ # if middle and lower are the same sign then lower bound = middle bound ==> change is between middle and upper
        a <- c # a is the lower bound value
        lower.fit <- mid.fit
        mid.init <- mid.fit$theta.hat
      }
      iter <- iter + 1
    }
    best.fit <- mid.fit
    return(best.fit)
  }


  if (h.method == '2step' & include.h==T){
    two.step.h <- T
  }

  # if initial values are not supplied, then they need to be randomly generated:
  if (is.null(init)){
    # if background risk is not included, then two.step.h is always false because starting value for background risk doesn't need to be generated
    if (include.h ==F) two.step.h <- F

    if (two.step.h & h.method=='2step'){
      # two.step.h = option to first generate starting values without background risk and then repeat with background risk
      if(!silent) print(noquote("Generating inital values without background risk."))
      # first we get initial values for parameters EXCEPT background risk (h)
      init.without.h <- init.generator(data, include.h =F, init.h.only =F, init=NULL, est.h=F, fixed.h = fixed.h)

      #init.without.h <- em.function.h(init.without.h, data, include.h=F)$theta.hat
      if (!silent) print(init.without.h)
      if(!silent) print(noquote("Generating inital values with background risk."))

      if(include.h) pars <- c('h', pars)
      init <- init.generator(data, include.h =T, init.h.only =T, init=init.without.h, fixed.h = fixed.h)
    }else if(h.method == 'bisection'){
      if(!silent) print(noquote("Running EM algorithm with bisection method for background risk."))
      upper.h <- 0.005

      lower.fit <- PICmodel.fit(l1_x , l2_x, pi_x, data, epsilon, short.epsilon, short.iter, short.runs, silent=T, init,
                              include.h=F, h.method="", two.step.h=F, include.priors, prior.type, fixed.h=NULL, intercept.prog, intercept.clear, intercept.prev)

      upper.fit <- PICmodel.fit(l1_x , l2_x, pi_x, data, epsilon, short.epsilon, short.iter, short.runs, silent=T,  init, #init = lower.fit$theta.hat,
                              include.h, h.method="", two.step.h=F, include.priors, prior.type, fixed.h=log(upper.h), intercept.prog, intercept.clear, intercept.prev)
      bisection.init <- suppressWarnings(bisection.em(0, upper.h, lower.fit, upper.fit)) # ignore warnings about negative sqrt
      if (! silent) print( c(h=bisection.init$fixed.h, bisection.init$theta.hat))
      final.res <- PICmodel.fit(l1_x , l2_x, pi_x, data, epsilon, short.epsilon, short.iter=1, short.runs=1, silent=F,  init= c(h=bisection.init$fixed.h, bisection.init$theta.hat),
                                include.h=T, h.method="", two.step.h=F, include.priors, prior.type, fixed.h=NULL, intercept.prog, intercept.clear, intercept.prev)
      return(final.res) # if bisection method just return results here
    }else{ # otherwise, generate initial values for all
      if(include.h & est.h) pars <- c('h', pars) # if background risk is included by initial values for h are not generated by 2-step then add h to pars
      init <- init.generator(data, include.h = include.h, init.h.only =F, init=NULL, est.h=est.h, fixed.h = fixed.h)
    }
  }else if(include.h & est.h) pars <- c('h', pars)


  # run the function
  if(!silent) print(noquote("Running EM algorithm."))
  final.res <- suppressWarnings(em.function.h(init, data, include.h, est.h, fixed.h))
  return(final.res)
}


