# Prevalence-incidence-cure model functions 
# 27/06/2024 added:
# - beta distribution for simulating screening time points
# - numerical differentiation for calculation of hessian after the EM algorithm has converged
# - fixed log-likelihood missing bracket mistake  
#### TO-DO: allow score test to evaluate a model with any covariates, currently it uses our specific 3 of covariates  

# model.fit(): function to fit our prevalence-incidence-cure model
# Inputs: 
# l1_x A vector containing the names of covariates used in the progression rate parameter (must match column name(s) in the input data), can be blank
# l2_x A vector containing the names of covariates used in the clearance rate parameter (must match column name(s) in the input data), can be blank
# pi_x A vector containing the names of covariates used in the parameter (probability of prevalent disease) (must match column name(s) in the input data), can be blank
# data Data used to fit the model containing columns for each term in l1_x, l2_x and pi_x. The first three columns must be (1)
# left interval, (2) right interval and (3) z indicator for prevalent/incident disease.The data must contain observed prevalent and incident cases, however some cases may be unknown (`NA`).
# num.runs Number of runs of the 'Short' EM algorithm used to determine initial values for the EM algorithm. Defaults to 30.
# short.epsilon Convergence criteria used in the 'Short' EM algorithm to determine initial values. Defaults to 0.1.
# epsilon Convergence criteria for the change in log-likelihood value used for stopping the EM algorithm. Defaults to 1e-08.
# silent Indicator variable for whether to print out each iteration or not. Default value set to FALSE.
# The output is a list containing the following elements:
#   initial.values - initial values determined by the Short EM algorithm process
#   theta.hat - optimum parameter values estimated by the EM algorithm
#   num.iterations - number of iterations until the EM algorithm converged
#   log.likelihood - value of the log.likelihood at the optimum parameter values
#   hess - hessian matrix
#   std.dev - standard deviation of parameter estimates
#   summary - data frame with estimate, std.dev, and 95\% CI for each parameter (useful for data set comparisons)
model.fit <- function(l1_x = c(), l2_x= c(), pi_x=c(), data, epsilon=1e-08, short.epsilon=1e-1, short.iter=10, short.runs=20,
                      silent=T,  init=NULL, include.h=T, two.step.h=T, include.priors=T, fixed.h=NULL, intercept.prog = T, intercept.clear = T, intercept.prev=T){
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
  log.likelihood.h <- function(current_par, data, include.h, est.h=T){
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
  
  estep.h <- function(current_par, data, include.h, est.h=T){
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
  
  build.expr.h <- function(model_par, include.h, est.h=T){
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
    
    p.expr <- ifelse(n3==1, str2expression(paste0("z*log(", expr3, ") + (1-z)*log(1-", expr3, ")")),
                     str2expression(paste0("z*(", expr3, ") - log(1 + exp(",expr3, "))")))
    f.expr <- str2expression(paste0("(1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",expr1 ,")/(exp(", expr1, ") + exp(",expr2,")))*(1-exp(-(exp(",expr1,") + exp(", expr2,"))*t))"))
    return(list(p.expr, f.expr))
  }
  
  mstep.h <- function(current_par, expected_z, data, include.h, est.h=T){
    z <- expected_z # expected value of z = output from the estep.h() function
    m1 <- n1+n2
    if (include.h){
      if (est.h) {
        m1 <- n1+n2 + 1
      }else {
        h <- fixed.h
      }
    }else {
      h <- -Inf
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
    exprs <- build.expr.h(pars, include.h, est.h)
    p.expr <- exprs[[1]] # prevalent function 
    f.expr <- exprs[[2]] # incident function
    
    # cauchy prior on the parameter: param = pars[i], mean=0, t=2.5 
    deriv.cauchy <- function(param, order=1){
      if(! include.priors) return(0)
      
      if (order==1){
        cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param, "^2 + 2.5^2))"))
        return(eval(D(cauchy.prior, param)))
      }else if (order==2){
        cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param[1], "^2 + 2.5^2))"))
        return(eval(D(D(cauchy.prior, param[1]), param[2])))
      }
    } 
    
    hess <- matrix(rep(0, (length(current_par))^2), nrow=length(current_par)) # empty hessian matrix
    grad <- rep(0, length(current_par)) # empty gradient matrix
    
    for (i in 1:m1){
      # we loop through the parameter list and calculate the first derivatives for gradient vector for the l1 and l2 parameters
      grad[i] <- sum(((1-z)*(deriv.f(right, pars[i]) - deriv.f(left, pars[i]))/(f(right) - f(left)))[z!=1]) + deriv.cauchy(pars[i])
      for (j in 1:i){
        # calculate 2nd derivatives for the Hessian matrix for h, l1 and l2 parameters; 
        # the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time 
        hess[i,j] <- hess[j,i] <- sum(((1-z)*(f(right) - f(left))^(-2)*
                                         ((deriv.f(right, c(pars[i], pars[j]), 2) - deriv.f(left, c(pars[i], pars[j]),2))*(f(right) - f(left)) -
                                            (deriv.f(right, pars[i], 1) - deriv.f(left, pars[i],1))*(deriv.f(right, pars[j],1) - deriv.f(left, pars[j], 1))))[z!=1]) + 
          deriv.cauchy(c(pars[i], pars[j]), order=2)
        
      }
    }
    # parameters for pi (for example p0, p1, p2, p3, p4) do not depend on left/right so we can calculate it separately
    for (i in (m1+1):length(current_par)){
      # calculate first derivatives for gradient vector for the pi parameters
      grad[i] <- sum(eval(D(p.expr, pars[i]))) + deriv.cauchy(pars[i])
      for (j in (m1+1):i){
        # calculate 2nd derivatives for Hessian matrix for the pi parameters; the hessian is symmetric so we can fill both hess[i,j] and hess[j,i] together to save time
        hess[i,j] <- hess[j,i] <- sum(eval(D(D(p.expr, pars[i]), pars[j]))) + deriv.cauchy(c(pars[i], pars[j]), order=2)
      }
    }
    new_par <-  current_par - solve(hess)%*%grad # single Newton step to calculate updated parameters
    return(list(as.vector(new_par), hess, grad)) # return the new parameters
  }
  
  # numerical.hess(): function to calculate the hessian using numerical methods
  numerical.hess <- function(mles, data, include.h, est.h=T){
    for (i in 1:length(pars)){    # assign the values of current_par to vectors with names from the 'pars' list (this will or won't have h and the number of parameters will correspond)
      assign(pars[i], mles[i], envir=environment())
    }
    
    deriv.cauchy <- function(param, order=1){
      if(! include.priors) return(0)
      if (order==1){
        cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param, "^2 + 2.5^2))"))
        return(eval(D(cauchy.prior, param)))
      }else if (order==2){
        cauchy.prior <- str2expression(paste0("log(2.5) - log(pi*(", param[1], "^2 + 2.5^2))"))
        return(eval(D(D(cauchy.prior, param[1]), param[2])))
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
        
        if (i ==j){ # only the diagonal entries have non-zero cauchy prior derivative
          hess[i,i] <- (-log.likelihood.h(mles+2*dh_i,  data, include.h, est.h) + 16*log.likelihood.h(mles+dh_i,  data, include.h, est.h) - 30*log.likelihood.h(mles,  data, include.h, est.h) +
                          16*log.likelihood.h(mles-dh_i,  data, include.h, est.h) - log.likelihood.h(mles-2*dh_i,  data, include.h, est.h))/(12*dh^2) + deriv.cauchy(c(pars[i], pars[i]), order=2)

        }else{ # accuracy to the order 4
          f_ii <- (-log.likelihood.h(mles+2*dh_i,  data, include.h, est.h) + 16*log.likelihood.h(mles+dh_i,  data, include.h, est.h) - 30*log.likelihood.h(mles,  data, include.h, est.h) +
                     16*log.likelihood.h(mles-dh_i,  data, include.h, est.h) - log.likelihood.h(mles-2*dh_i,  data, include.h, est.h))/(12*dh^2) + deriv.cauchy(c(pars[i], pars[i]), order=2)
          f_jj <- (-log.likelihood.h(mles+2*dh_j,  data, include.h, est.h) + 16*log.likelihood.h(mles+dh_j,  data, include.h, est.h) - 30*log.likelihood.h(mles,  data, include.h, est.h) +
                     16*log.likelihood.h(mles-dh_j,  data, include.h, est.h) - log.likelihood.h(mles-2*dh_j,  data, include.h, est.h))/(12*dh^2) + deriv.cauchy(c(pars[j], pars[j]), order=2)

          hess_entry <- (16*(log.likelihood.h(mles+dh_ij,  data, include.h, est.h) + log.likelihood.h(mles-dh_ij,  data, include.h, est.h)) -
                           (log.likelihood.h(mles+2*dh_ij,  data, include.h, est.h) + log.likelihood.h(mles-2*dh_ij,  data, include.h, est.h)) -
                           30*log.likelihood.h(mles,  data, include.h, est.h))/(24*dh^2) - f_ii/2 - f_jj/2

          hess[i, j] <- hess[j, i] <- hess_entry
        }
      }
    }
    return(hess)
  }
  
  # em.function.h(): combining the E- and M-step and repeating until the parameter estimates converge
  em.function.h <- function(init, data, include.h, est.h=T){
    new_theta <- init
    old_llk <- 0
    new_llk <- 100 # make big to start with
    iter <- 1
    
    while ((abs(new_llk - old_llk) > epsilon)){
      current_theta <- new_theta
      old_llk <- new_llk
      next_em_step <- mstep.h(current_theta, estep.h(current_theta, data, include.h, est.h), data, include.h, est.h)
      new_theta <- as.vector(next_em_step[[1]])
      new_llk <- log.likelihood.h(new_theta, data, include.h, est.h)
      if (any(abs(current_theta - new_theta) > 1)){
        print(current_theta)
        
      } 
      if (!silent & (iter %/% 10 ==0 | iter %% 10 ==0)) print(round(c(new_theta, new_llk), 4))
      iter <- iter + 1
      if (iter > 1000){
        break
      }
    }
    hess <- numerical.hess(new_theta, data, include.h, est.h) 
    names(new_theta) <- pars
    std.dev <- round(sqrt(diag(solve(-hess))),4)
    param <- c(if(include.h){"h"}, if(intercept.prog) {"intercept"}, l1_x, if(intercept.clear) {"intercept"}, l2_x, if(intercept.prev) {"intercept"}, pi_x)
    summary <- round(data.frame(theta.hat = new_theta, std.dev, lower = new_theta - 1.96*std.dev, upper = new_theta + 1.96*std.dev), 4)
    summary <- cbind(param, summary)
    rownames(summary) <- names(new_theta)
    return(list(model = list(l1_x, l2_x, pi_x), initial.values = init, fixed.h = fixed.h, theta.hat = new_theta, num.iterations=iter, 
                log.likelihood = new_llk, hess=hess, grad = next_em_step[[3]], std.dev=std.dev, summary=summary))
  }
  
  # short.em():  lots of short runs of the em function to generate appropriate starting values
  short.em <- function(data, init.h.only, include.h, init, est.h=T){
    if(!init.h.only){ #initial values for only l1, l2, p
      new_theta <- log(runif(n1+n2+n3))
      if (n3==1 & intercept.prev) new_theta[n1+n2+n3] <- exp(new_theta[n1+n2+n3]) # if only one parameter for prevalence, then no covariates then they're not on log scale  
      if(include.h & est.h) { # only generate if estimating h 
        new_theta <- c(log(runif(1, 0, 0.02)), new_theta)
      }
    }else{# initial values for background risk using previously found values of l1, l2, p
      new_theta <- c(log(runif(1, 0, 0.02)), init)
    }
    new_llk <- 100
    old_llk <- 0
    iter <- 1
    while ((abs(new_llk - old_llk) > short.epsilon) & iter < short.iter){
      current_theta <- new_theta
      old_llk <- new_llk
      iter <- iter + 1
      new_theta <- try(as.vector(mstep.h(current_theta, estep.h(current_theta, data, include.h, est.h), data, include.h, est.h)[[1]]), silent=T)
      if (class(new_theta) == 'try-error') return(short.em(data, init.h.only, include.h, init, est.h))
      new_llk <- log.likelihood.h(new_theta, data, include.h, est.h)
      if (abs(new_llk) == Inf | is.nan(new_llk) | is.na(new_llk)) return(short.em(data, init.h.only, include.h, init, est.h))
    }
    return(list(theta.hat = new_theta, log.likelihood = new_llk))
  }
  
  # init.generator(): performs short runs of the em function using short.em and then returns the most common set of 
  # starting values with the highest loglikelihood as the initial values
  init.generator <- function(data, init.h.only, include.h, init, est.h=T){
    if (!silent) pb <- txtProgressBar(min = 0, max = short.runs, style = 3, width = 50, char = "=")
    short.inits <- list()
    for(k in 1:short.runs) {
      short.inits[[k]] <- short.em(data, init.h.only, include.h, init, est.h)[c("theta.hat", "log.likelihood")]
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
  
  if (!is.null(fixed.h)) {
    two.step.h <- F
    include.h <- T
    est.h <- F
  }else{
    est.h <- T
  }
  
  # if initial values are not supplied, then they need to be randomly generated:
  if (is.null(init)){
    # if background risk is not included, then two.step.h is always false because starting value for background risk doesn't need to be generated 
    if (include.h ==F) two.step.h <- F
    
    if (two.step.h){
      # two.step.h = option to first generate starting values without background risk and then repeat with background risk 
      if(!silent) print(noquote("Generating inital values without background risk."))
      # first we get initial values for parameters EXCEPT background risk (h)
      init.without.h <- init.generator(data, include.h =F, init.h.only =F, init=NULL)
      
      #init.without.h <- em.function.h(init.without.h, data, include.h=F)$theta.hat
      if (!silent) print(init.without.h)
      if(!silent) print(noquote("Generating inital values with background risk."))
      
      if(include.h) pars <- c('h', pars) 
      init <- init.generator(data, include.h =T, init.h.only =T, init=init.without.h)
      
    }else{ # otherwise, generate initial values
      if(include.h & est.h) pars <- c('h', pars) # if background risk is included by initial values for h are not generated by 2-step then add h to pars
      init <- init.generator(data, include.h = include.h, init.h.only =F, init=NULL, est.h=est.h)
    }
  }
  
  # run the function
  if(!silent) print(noquote("Running EM algorithm."))
  final.res <- em.function.h(init, data, include.h, est.h)
  return(final.res)
}



# gamma.simulator(): function to simulate data from our model for a specified shape parameter (k), screening interval, with/without background risk (include.h)
# note: 
#   - parameters go in order (background risk, lambda_1, lambda_2, prevalent)
#   - this code is simplified to not include covariates    
gamma.simulator <- function(n, params,  k, interval=3, include.h=T){
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



library(stringr)
# gamma.score.test(): function to compute p-value of score test given model fit, step size (dh), order (accuracy, 1, 2 or 4) and data set
gamma.score.test <- function(fit, dh, accuracy, data){
  require("stringr")
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
