% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_predict.R
\name{PICmodel.predict}
\alias{PICmodel.predict}
\title{Make predictions from the prevalence-incidence-cure model}
\usage{
PICmodel.predict(
  l1_x,
  l2_x,
  pi_x,
  data,
  time.points,
  fit,
  calc.CI = F,
  include.h = T
)
}
\arguments{
\item{l1_x}{A vector containing the names of covariates used in the \ifelse{html}{\out{\eqn{\lambda}<sub>1</sub>}}{ \eqn{\lambda_1}} (progression rate) parameter (must match column name(s) in the input data)}

\item{l2_x}{A vector containing the names of covariates used in the \eqn{\lambda_2} (clearance rate) parameter (must match column name(s) in the input data)}

\item{pi_x}{A vector containing the names of covariates used in the \eqn{\pi} parameter (probability of prevalent disease) (must match column name(s) in the input data)}

\item{data}{Data set of covariates from which to make the predictions.}

\item{time.points}{Numeric vector of time points used to make cumulative risk predictions}

\item{fit}{Parameter estimates for the model to be used (output from model.fit)}

\item{calc.CI}{Indicator variable for whether confidence intervals for the cumulative risk should be calculated. Defaults to FALSE.}

\item{include.h}{Indicator variable for whether background risk was included in the model fitting procedure. Defaults to TRUE}
}
\value{
The output is a list for each unique combination of covariates used in the model containing the time points, the cumulative risk estimate along with the upper and lower 95% confidence intervals.
}
\description{
Make predictions from the prevalence-incidence-cure model
}
\author{
Kelsi Kroon, Hans Bogaards, Hans Berkhof
}
