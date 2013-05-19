\name{SODModel}
\alias{SODModel}
\title{Run the SOD model}
\usage{
  SODModel(parms, times, locations, init, reps = 1,
    lambda.ex = 0, K = 50, stochastic.e = NULL,
    stochastic.d = FALSE, verbose = interactive(),
    parallel = FALSE)
}
\arguments{
  \item{parallel}{Run the simulation in parallel or not?}

  \item{...}{arguments to be passed to SODModel}
}
\description{
  Run the SOD model
}

