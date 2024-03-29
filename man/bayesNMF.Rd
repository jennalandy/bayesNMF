% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesNMF.R
\name{bayesNMF}
\alias{bayesNMF}
\title{Perform single-study Bayesian NMF with the provided likelihood and prior
combination. Exact rank `N` or maximum rank `max_N` must be provided.}
\usage{
bayesNMF(
  M,
  N = NULL,
  max_N = NULL,
  inits = NULL,
  likelihood = "normal",
  prior = "truncnormal",
  prior_parameters = NULL,
  logevery = 100,
  file = paste0("nmf_", likelihood, "_", prior),
  overwrite = FALSE,
  true_P = NULL,
  niters = NULL,
  burn_in = NULL
)
}
\arguments{
\item{M}{mutational catalog matrix, K x G}

\item{N}{number of signatures}

\item{prior}{string, one of c('truncnormal','exponential')}

\item{logevery}{the log, save, and plot files will be updated every
`logevery` iterations}

\item{file}{file name without extension of log, save, and plot files}

\item{overwrite}{if `overwrite = TRUE`, the log, safe, and plot files of
previous runs with the same `file` will be overwritten}

\item{true_P}{(optional) true signatures matrix P to compare to in a heatmap}

\item{niters}{how many iterations the Gibbs sampler is run}

\item{burn_in}{the first `burn_in` iterations will be ignored when computing
MAP estimate}

\item{P}{(optional) initial signatures matrix, K x N}

\item{E}{(optional) initial exposure matrix, N x G}

\item{sigmasq}{(optinal) initial variance vector, length K}
}
\value{
list
}
\description{
Perform single-study Bayesian NMF with the provided likelihood and prior
combination. Exact rank `N` or maximum rank `max_N` must be provided.
}
