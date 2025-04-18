% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_gamma.R
\name{score.test.gamma}
\alias{score.test.gamma}
\title{Gamma Score test}
\usage{
score.test.gamma(fit, dh = 1e-05, accuracy = 4, data)
}
\arguments{
\item{fit}{Output from model fit.}

\item{dh}{Step size used in numerical differentiation (e.g., 1e-5, 1e-10)}

\item{accuracy}{Order of accuracy in numerical differentiation, options are 1, 2, or 4.}

\item{data}{Data used to fit the model.}
}
\value{
The output is a list containing the following elements:
\itemize{
\item grad: value of the gradient at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
\item inv.hess: value of the inverse hessian at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
\item hess: value of the hessian at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
\item loglik: value of the log-likelihood at the shape parameter equal to 1 and the rest of the parameters at their maximum likelihood estimates.
\item p.val: p-value of the score test for whether Cause 1 follows a Gamma rather than Exponential distribution (tests whether shape parameter equals 1 under \eqn{H_0} or not).
}
}
\description{
Function to compute p-value of the score test given model fit, step size (dh), order (accuracy, 1, 2 or 4) and data.
}
\author{
Kelsi Kroon, Hans Bogaards, Hans Berkhof
}
