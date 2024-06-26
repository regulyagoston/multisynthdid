% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver.R
\name{sc.weight.fw.covariates}
\alias{sc.weight.fw.covariates}
\title{A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0}
\usage{
sc.weight.fw.covariates(
  Y,
  theta,
  X,
  beta_valid,
  zeta.lambda = NULL,
  zeta.omega = NULL,
  lambda.intercept = TRUE,
  omega.intercept = TRUE,
  min.decrease = 0.001,
  max.iter = 1000,
  lambda = NULL,
  omega = NULL,
  beta = NULL,
  update.lambda = TRUE,
  update.omega = TRUE
)
}
\description{
A Frank-Wolfe + Gradient solver for lambda, omega, and beta when there are covariates
Uses the exact line search Frank-Wolfe steps for lambda, omega and (1/t)*gradient steps for beta
pass update.lambda=FALSE/update.omega=FALSE to fix those weights at initial values, defaulting to uniform 1/T0 and 1/N0
}
