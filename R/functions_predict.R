#' Make predictions from the prevalence-incidence-cure model
#' @param l1_x A vector containing the names of covariates used in the \ifelse{html}{\out{\eqn{\lambda}<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter (must match column name(s) in the input data)
#' @param l2_x A vector containing the names of covariates used in the \eqn{\lambda_2} (clearance rate) parameter (must match column name(s) in the input data)
#' @param pi_x A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease) (must match column name(s) in the input data)
#' @param data Data set of covariates from which to make the predictions.
#' @param time.points Numeric vector of time points used to make cumulative risk predictions
#' @param fit Parameter estimates for the model to be used (output from model.fit)
#' @param calc.CI Indicator variable for whether confidence intervals for the cumulative risk should be calculated. Defaults to FALSE.
#' @param include.h Indicator variable for whether background risk was included in the model fitting procedure. Defaults to TRUE
#' @return The output is a list for each unique combination of covariates used in the model containing the time points, the cumulative risk estimate along with the upper and lower 95% confidence intervals.
#'
#' @author Kelsi Kroon, Hans Bogaards, Hans Berkhof
#' @export
model.predict <- function(l1_x, l2_x, pi_x, data, time.points, fit, calc.CI=F, include.h=T){
  g <- function(l1, l2, t){
    return(l1 / (l1 + l2) * (1 - exp(-(l1 + l2)* t)))
  }

  se.cumrisk <- function(t, data1, data2, data3, n1, n2, n3, fit, include.h){
    create.par <- function(n1, n2, n3) {
      par1 <- paste0(c("g"), seq(0, n1 - 1))
      par2 <- paste0(c("w"), seq(0, n2 - 1))
      par3 <- paste0(c("p"), seq(0, n3 - 1))
      return(c(par1, par2, par3))
    }

    if (include.h){
      pars <- c("h", create.par(n1, n2, n3))
    }else{
      pars <- create.par(n1, n2, n3)
      h <- 0
    }

    if (n1 !=1){
      for (i in 1:(n1-1)){
        assign(paste0("data1", i), data1[,i+1], envir=environment())
      }}
    if (n2 !=1){
      for (i in 1:(n2-1)){
        assign(paste0("data2", i), data2[,i+1], envir=environment())
      }}
    if (n3!=1){
      for (i in 1:(n3-1)){
        assign(paste0("data3", i), data3[,i+1], envir=environment())
      }}


    # assign the values of current_par to vectors with names from the 'pars' list
    for (i in 1:length(pars)){
      assign(pars[i], fit[["theta.hat"]][i], envir=environment())
    }


    # create the expression for the complementary log log cumulative risk, it is coded this way to allow for user-defined covariates
    # and then we take the derivative so it needs to be in expression form
    create.cllcr.expr <- function(n1, n2, n3, pars, data1, data2, data3){
      pars <- pars[-1]
      # expression for lambda_1 (progression rate parameter)
      expr1 <- pars[1]
      if (n1 != 1){
        for (i in 1:(n1-1)){
          expr1 <- paste0(expr1, "+", pars[1+i], "*", "data1", i)
        }}

      # expression for lambda_2 (clearance rate parameter)
      expr2 <- pars[n1+1]
      if (n2!=1){
        for (i in 1:(n2-1)){
          expr2 <- paste0(expr2, "+", pars[n1+1+i], "*", "data2", i)
        }}

      # expression for pi (prevalent probability parameter)
      expr3 <- pars[n1+n2+1]
      if (n3 > 1){
        for (i in 1:(n3-1)){
          expr3 <- paste0(expr3, "+", pars[n1+n2+1+i], "*", "data3", i)
        }

        #  expected complete log-likelihood for Ri!=Inf
        llcr <- paste0("log(- log( (exp(", expr3, ")/(1+ exp(", expr3, "))) + (1/(1+exp(", expr3, ")))*((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
                       expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t)))))" )


      }else if (n3 == 1){ # if there are no covariates for pi we use a different expression (don't need to use logit function)
        #  expected complete log-likelihood for Ri!=Inf for when there is only intercept for pi parameter
        # llcr <- paste0("log(- log( log(", expr3, ") + (log(1-", expr3, ") + log((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
        #                expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t))))))" )
        llcr <- paste0("log(- log( ", expr3, " + (1 - ", expr3, ")*((1-exp(-exp(h)*t)) + exp(-exp(h)*t)*(exp(",
                       expr1, ")/(exp(", expr1, ") + exp(", expr2, ")))*(1-exp(-(exp(", expr1, ") + exp(", expr2,"))*t)))))" )
      }

      # convert from string to expression to be used in the differentiate function in the gradient calculation
      return(str2expression(llcr))
    }

    # function to calculate the gradient of the log log cumulative risk using the built in R derivative function D()
    grad.loglog <- function(par){
      return(eval(D(create.cllcr.expr(n1, n2, n3, pars, data1, data2, data3), par)))
    }

    return(sqrt(as.numeric(t(sapply(pars, grad.loglog)) %*% (solve(-fit[["hess"]])) %*% (sapply(pars, grad.loglog)))))
  }


  requireNamespace("dlpyr", quietly = TRUE)
  theta.hat <- fit$theta.hat
  if (include.h) {
    h <- exp(theta.hat[1])
    theta.hat <- theta.hat[-1]
  }else{
    h <- 0
  }

  n1 <- length(l1_x) + 1
  n2 <- length(l2_x) + 1
  n3 <- length(pi_x) + 1

  preds <- list()
  data <- unique(data)
  rownames(data) <- seq(1:nrow(data))
  # calculate the cumulative risk for each unique combination of covariates in the input prediction data and store in a list
  for (row in 1:nrow(data)){
    sub.data <- subset(data,  subset=rownames(data) == row)

    create.covariate.data <- function(data){
      data[['intercept']] <- rep(1, length(dim(data)[1]))
      data1 <- as.matrix(data[,c("intercept", l1_x)])
      data2 <- as.matrix(data[,c("intercept", l2_x)])
      data3 <- as.matrix(data[,c("intercept", pi_x)])
      return(list(data1, data2, data3))
    }

    covariate.data <- create.covariate.data(sub.data)
    data1 <- covariate.data[[1]]
    data2 <- covariate.data[[2]]
    data3 <- covariate.data[[3]]

    l1 <- as.numeric(exp(data1 %*% theta.hat[1:(n1)]))
    l2 <- as.numeric(exp(data2 %*% theta.hat[(n1+1):(n1+n2)]))
    if (n3 == 1){
      #if pi has no covariates then the value of p is the same for all women
      p <- theta.hat[n1+n2+n3]
    }else if (n3 > 1){
      p <- as.numeric(exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])/(1+exp(data3 %*% theta.hat[(n1+n2+1):(n1+n2+n3)])))
    }
    prob <- p + (1-p)*(1 - exp(-h*time.points) + exp(-h*time.points)*g(l1, l2, time.points))

    ## calculate confidence errors for cumulative risk using the complementary log-log cumulative risk
    if (calc.CI){
      cr.est <- log(-log(prob))
      cr.se <- sapply(time.points,  se.cumrisk, data1, data2, data3, n1, n2, n3, fit, include.h)
      cr.lci <- exp(-exp(cr.est + 1.96*cr.se))
      cr.uci <- exp(-exp(cr.est - 1.96*cr.se))
      preds[[row]] <- data.frame(Time=time.points, CR = prob, CR.se = cr.se, CR.lower95 = cr.lci, CR.upper95 = cr.uci)
    }else{
      preds[[row]] <- data.frame(Time=time.points, CR = prob)
    }

  }
  return(preds)
}
