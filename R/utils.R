# import pipe operator from magrittr

#' Pipe operator
#'
#' See \code{\link[magrittr:pipe]{%>%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
NULL

# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PUBLIC METHODS ##############
####################################


#' Compute the expected data matrix
#' @param self bayesNMF_sampler object
#' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
#' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
#' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
#' 
#' @return matrix, dimensions K x G
#' @noRd
get_Mhat_ <- function(self, P = NULL, A = NULL, E = NULL) {
  # use current parameters if not provided
  if (is.null(P)) { P <- self$params$P }
  if (is.null(A)) { A <- self$params$A }
  if (is.null(E)) { E <- self$params$E }

  # handle single signature case
  if (self$dims$N == 1) {
    # convert P to a column vector
    P <- matrix(P, nrow = self$dims$K, ncol = 1)
    # Keep A as a 1x1 matrix, not a scalar
    diag_A <- A
  } else {
    # convert A to a diagonal matrix
    diag_A <- diag(A[1, ])
  }

  # compute matrix product P*A*E
  Mhat <- P %*% diag_A %*% E
  return(Mhat)
}

#' Compute the log likelihood
#' @param self bayesNMF_sampler object
#' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
#' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
#' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
#' @param sigmasq vector, optional, dimensions G, uses self$params$sigmasq if not provided, only used if likelihood is "normal"
#' @param likelihood string, likelihood distribution for M, one of "poisson", "normal", uses self$specs$likelihood if not provided, created as an option to compute acceptance ratios in MH updates
#' @param return_matrix boolean, whether to return the log likelihood matrix, used to compute acceptance ratios in MH updates (default FALSE)
#' 
#' @return scalar log likelihood if return_matrix is FALSE, otherwise matrix of log likelihoods per observed value (K x G)
#' @noRd
get_loglik_ <- function(self, P = NULL, A = NULL, E = NULL, sigmasq = NULL, likelihood = self$specs$likelihood, return_matrix = FALSE) {
  # use current parameters in self if not provided
  if (is.null(P)) { P <- self$params$P }
  if (is.null(A)) { A <- self$params$A }
  if (is.null(E)) { E <- self$params$E }

  # compute expected data matrix given parameters
  Mhat <- self$get_Mhat(P = P, A = A, E = E)

  # compute log likelihood
  if (likelihood == "normal") {
    if (is.null(sigmasq)) {
      sigmasq <- self$params$sigmasq
    }

    # if a variance of length G is provided, convert to matrix
    # with each column being the same value sigmasq_g
    if (!("matrix" %in% class(sigmasq))) {
      sigmasq <- matrix(
        rep(sigmasq, self$dims$K),
        ncol = self$dims$G,
        nrow = self$dims$K,
        byrow = TRUE
      )
    }

    # dnorm doesn't like matrices, so we need to loop over columns
    loglik_matrix <- lapply(1:self$dims$G, function(g) {
      dnorm(
        self$data[, g],
        mean = Mhat[, g],
        sd = sqrt(sigmasq[, g]),
        log = TRUE
      ) # K x 1 column vector
    }) %>%
      do.call(cbind, .) # K x G matrix
  } else if (likelihood == "poisson") {
    # clip Mhat to avoid log(0)
    Mhat <- pmax(Mhat, 1e-6)
    # loop over columns of data and Mhat to compute log likelihood
    loglik_matrix <- lapply(1:self$dims$G, function(g) {
      dpois(self$data[, g], lambda = Mhat[, g], log = TRUE) # K x 1 column vector
    }) %>%
      do.call(cbind, .) # K x G matrix
  }
  if (return_matrix) {
    return(loglik_matrix)
  } else {
    return(sum(loglik_matrix))
  }
}

#' Compute the log posterior using model specifications
#' @param self bayesNMF_sampler object
#' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
#' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
#' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
#' @param sigmasq vector, optional, dimensions G, uses self$params$sigmasq if not provided, only used if likelihood is "normal"
#' 
#' @return scalar log posterior
#' @noRd
get_logpost_ <- function(self, P = NULL, A = NULL, E = NULL, sigmasq = NULL) {
  # use current parameters in self if not provided
  if (is.null(P)) { P <- self$params$P }
  if (is.null(A)) { A <- self$params$A }
  if (is.null(E)) { E <- self$params$E }
  if (is.null(sigmasq)) { sigmasq <- self$params$sigmasq }

  # compute log prior
  logprior <- 0
  if (self$specs$prior == "truncnormal") {
    for (n in 1:self$dims$N) {
      logprior <- logprior + truncnorm::dtruncnorm(
        P[, n],
        mean = self$prior_params$Mu_p[, n],
        sd = sqrt(self$prior_params$Sigmasq_p[, n]),
        a = 0, b = Inf
      ) %>% log() %>% sum()
      logprior <- logprior + truncnorm::dtruncnorm(
        E[n, ],
        mean = self$prior_params$Mu_e[n, ],
        sd = sqrt(self$prior_params$Sigmasq_e[n, ]),
        a = 0, b = Inf
      ) %>% log() %>% sum()
    }
  } else if (self$specs$prior == "exponential") {
    for (n in 1:self$dims$N) {
      logprior <- logprior + dexp(
        P[, n],
        rate = self$prior_params$Lambda_p[, n],
        log = TRUE
      ) %>% sum()
      logprior <- logprior + dexp(
        E[n, ],
        rate = self$prior_params$Lambda_e[n, ],
        log = TRUE
      ) %>% sum()
    }
  } else if (self$specs$prior == "gamma") {
    for (n in 1:self$dims$N) {
      logprior <- logprior + dgamma(
        P[, n],
        shape = self$prior_params$Alpha_p[, n],
        rate = self$prior_params$Beta_p[, n],
        log = TRUE
      ) %>% sum()
      logprior <- logprior + dgamma(
        E[n, ],
        shape = self$prior_params$Alpha_e[n, ],
        rate = self$prior_params$Beta_e[n, ],
        log = TRUE
      ) %>% sum()
    }
  }

  # compute log likelihood
  # defaults to specified likelihood in self$specs$likelihood
  loglik <- self$get_loglik(P, A, E, sigmasq)

  # return log posterior
  return(loglik + logprior)
}

#' Compute the maximum a posteriori (MAP) estimate
#' @param self bayesNMF_sampler object
#' @param end_iter integer, last iteration to consider for inference, defaults to current iteration
#' @param n_samples integer, number of samples to consider for inference, defaults to MAP_over specified in convergence control
#' @param final boolean, if TRUE, subset to only included signatures
#' @param credible_interval float, credible interval width (default 0.95)
#' 
#' @return None, updates self$MAP and self$credible_intervals
#' @noRd
get_MAP_ <- function(
  self,
  end_iter = self$state$iter,
  n_samples = self$specs$convergence_control$MAP_over,
  final = FALSE,
  credible_interval = 0.95
) {
  # across the past n_samples samples,
  # i. find mode of A matrix
  # among samples that match A mode
  #   ii. renormalize so that columns of P sum to 1
  #   iii. compute element-wise average of P and E
  # iv. return MAP P, A, E

  # define indices to consider
  # end_iter can only be provided if self$specs$save_all_samples
  if (!self$specs$save_all_samples & end_iter != self$state$iter) {
    stop("end_iter cannot be provided unless self$specs$save_all_samples is TRUE")
  }
  if (end_iter == self$state$iter) {
    considered_idx <- self$state$MAP_idx
  } else {
    considered_idx <- seq(
      from = end_iter - n_samples + 1,
      to = end_iter
    )
  }

  # i. find mode of A matrix
  A_MAP = get_mode(self$samples$A[considered_idx])
  MAP_idx = considered_idx[A_MAP$idx]
  if (final) {
    keep_sigs <- which(A_MAP$matrix[1,] == 1)
  } else {
    keep_sigs <- 1:ncol(A_MAP$matrix)
  }
  A_MAP$matrix <- A_MAP$matrix[, keep_sigs, drop = FALSE]

  #   ii. renormalize so that columns of P sum to 1
  renormalized_samples <- lapply(MAP_idx, function(idx) {
    renormalized <- renormalize(
      self$samples$P[[idx]], self$samples$E[[idx]]
    )

    return(list(
      P = renormalized$P[, keep_sigs, drop = FALSE],
      E = renormalized$E[keep_sigs, , drop = FALSE]
    ))
  })

  #   iii. compute element-wise average of P and E
  P_MAP <- Reduce("+", lapply(renormalized_samples, `[[`, "P")) / 
    length(renormalized_samples)
  E_MAP <- Reduce("+", lapply(renormalized_samples, `[[`, "E")) / 
    length(renormalized_samples)

  # iv. return MAP P, A, E
  self$MAP <- list(
    P = P_MAP,
    A = A_MAP$matrix,
    E = E_MAP,
    idx = MAP_idx,
    A_counts = A_MAP$top_counts,
    keep_sigs = keep_sigs
  )

  # if final, compute credible interval for P and E
  # or to update credible intervals if they have been previously computed
  if (final | !is.null(self$credible_intervals)) {
    probs <- c(0.5 - credible_interval / 2, 0.5 + credible_interval / 2)
    P_arr <- abind::abind(lapply(renormalized_samples, `[[`, "P"), along = 3)
    P_CI <- apply(P_arr, c(1, 2), quantile, probs = probs)
    P_CI_list <- list(
      lower = P_CI[1, , ],
      upper = P_CI[2, , ]
    )
    E_arr <- abind::abind(lapply(renormalized_samples, `[[`, "E"), along = 3)
    E_CI <- apply(E_arr, c(1, 2), quantile, probs = probs)
    E_CI_list <- list(
      lower = E_CI[1, , ],
      upper = E_CI[2, , ]
    )
    self$credible_intervals <- list(
      P = P_CI_list,
      E = E_CI_list
    )
  }
}




####################################
###### PRIVATE METHODS #############
####################################


#' Get tempering schedule
#' @description
#' Temperature parameter is used to control the balance between exploration and exploitation in the sampler.
#' The temperature starts at 0 and slowly increases to 1 over the course of the first `n_temp` iterations.
#' 
#' @param len integer, total number of iterations
#' @param n_temp integer, number of temperature levels
#' @return vector of temperatures, length `len`
#' @noRd
get_temp_sched_ <- function(len, n_temp) {
  nX <- round(n_temp / 374) # number of iterations per temperature
  nX <- max(nX, 1)
  temp_sched <- c(
    rep(0, nX), # initial temperature 0
    c(sapply(9:5, function(x) {rep(10^(-x), nX)})),
    rep(10**(-4), round(8 * nX)),
    c(sapply(4:1, function(y) {
      c(sapply(seq(0, 8.9, by = 0.1), function(x) {
        rep((1+x)*10**(-y), nX)
      }))
    }))
  )

  # if already too long, sample n_temp in order
  if (length(temp_sched) > n_temp) {
    temp_sched <- sort(sample(temp_sched, n_temp))
  }

  # add final temperature (1) to the end
  temp_sched <- c(
    temp_sched,
    rep(1, len - length(temp_sched)) # final temperature at 1
  )
  return(temp_sched)
}

#' Update sample metrics
#' @param self bayesNMF_sampler object
#' @param update_trace boolean, whether to update the trace plot (default FALSE)
#' @return None, updates self$state$sample_metrics and saves trace plot if update_trace is TRUE
#' @noRd
update_sample_metrics_ <- function(self, update_trace = FALSE) {
  metrics_i <- compute_metrics_(self, final = FALSE)
  metrics_i$rank <- sum(self$params$A)
  metrics_i$temp <- self$temperature_schedule[self$state$iter]
  self$state$sample_metrics <- rbind(self$state$sample_metrics, metrics_i)
  if (update_trace) {
    options(bitmapType = "cairo")
    trace_plot(self, save = TRUE)
  }
}

#' Update MAP metrics
#' @param self bayesNMF_sampler object
#' @param final boolean, if TRUE, subset to only included signatures
#' 
#' @return None, updates self$state$MAP_metrics and saves trace plot if periodic save or final is TRUE
#' @noRd
update_MAP_metrics_ <- function(self, final = FALSE) {
  metrics_i <- compute_metrics_(
    self, P = self$MAP$P, A = self$MAP$A, E = self$MAP$E, final = final, MAP = TRUE
  )
  options(bitmapType = "cairo")

  MAP_idx <- self$state$MAP_idx
  if (!self$specs$save_all_samples) {
    # if not saving all samples, MAP_idx is between 1 and MAP_over
    # need to scale up so that MAP_over is moved to current iteration
    MAP_idx <- MAP_idx + self$state$iter - self$specs$convergence_control$MAP_over
  } 

  # for MAP metrics, likelihood and posterior need to be computed per-MAP sample then averaged
  # because P and E have been renormalized so prior no longer applies
  # subset sample_metrics to MAP idx, then average
  metrics_i$loglikelihood <- self$state$sample_metrics %>%
    dplyr::filter(iter %in% MAP_idx) %>%
    dplyr::pull(loglikelihood) %>%
    mean()
  metrics_i$logposterior <- self$state$sample_metrics %>%
    dplyr::filter(iter %in% MAP_idx) %>%
    dplyr::pull(logposterior) %>%
    mean()

  # recompute BIC with this likelihood
  metrics_i$BIC <- -2 * metrics_i$loglikelihood +
    metrics_i$n_params * log(self$dims$G)

  # add rank and MAP_A counts
  metrics_i$rank <- sum(self$MAP$A)
  metrics_i$MAP_A_counts <- unname(self$MAP$A_counts[1])
  metrics_i$mean_temp <- mean(self$temperature_schedule[MAP_idx])
  
  # update MAP metrics
  self$state$MAP_metrics <- rbind(self$state$MAP_metrics, metrics_i)

  # save trace plot if periodic save or final
  if (self$specs$periodic_save | final) { 
    trace_plot(self, MAP_means = TRUE, save = TRUE)
  }
}

#' Compute metrics for a sample
#' @description
#' Helper for update_sample_metrics_ and update_MAP_metrics_
#' 
#' @param self bayesNMF_sampler object
#' @param P matrix, optional, dimensions K x N, uses self$params$P if not provided
#' @param A matrix, optional, dimensions 1 x N, uses self$params$A if not provided
#' @param E matrix, optional, dimensions N x G, uses self$params$E if not provided
#' @param final boolean, if TRUE, subset to only included signatures
#' @param MAP boolean, if TRUE, ignores loglikelihood and logposterior calculations because P and E have been renormalized so prior no longer applies
#' 
#' @return list of metrics, including iter, RMSE, KL, loglikelihood, logposterior, n_params, and BIC
#' @noRd
compute_metrics_ <- function(self, P = NULL, A = NULL, E = NULL, final = FALSE, MAP = FALSE) {
  if (is.null(P)) { P <- self$params$P }
  if (is.null(A)) { A <- self$params$A }
  if (is.null(E)) { E <- self$params$E }
  
  # if final, P and E have already been filtered to included signatures
  # recode A as a diagonal of all 1s
  if (final) {
    A <- matrix(1, nrow = 1, ncol = ncol(P))
  }

  Mhat <- self$get_Mhat(P = P, A = A, E = E)
  n_params <- sum(A) * (self$dims$G + self$dims$K)
  if (MAP) {
    loglikelihood <- NA_real_
    logposterior <- NA_real_
    BIC <- NA_real_
  } else {
    loglikelihood <- self$get_loglik(P = P, A = A, E = E)
    logposterior <- self$get_logpost(P = P, A = A, E = E)
    BIC <- -2 * loglikelihood + n_params * log(self$dims$G)
  }
  
  metrics <- list(
    iter = self$state$iter,
    RMSE = sqrt(mean((Mhat - self$data) ** 2)),
    KL = padded_KL_(Mhat, self$data),
    loglikelihood = loglikelihood,
    logposterior = logposterior,
    n_params = n_params,
    BIC = BIC
  )
  if (self$specs$MH) {
    # average acceptance rates for active signatures
    metrics$P_mean_acceptance_rate <- mean(
      self$acceptance_rates$P_acceptance_rate[, self$params$A[1,] == 1]
    )
    metrics$E_mean_acceptance_rate <- mean(
      self$acceptance_rates$E_acceptance_rate[self$params$A[1,] == 1, ]
    )
  }

  return(metrics)
}

#' Compute KL divergence with padding to avoid log(0)
#' @description
#' Helper for compute_metrics_
#' Handles cases where Mhat or M is 0 by adding a small amount to both.
#' 
#' @param Mhat matrix, dimensions K x G
#' @param M matrix, dimensions K x G
#' 
#' @return scalar KL divergence
#' @noRd
padded_KL_ <- function(Mhat, M) {
  Mhat <- pmax(Mhat, 1e-6)
  M <- pmax(M, 1e-6)
  return(sum(M * log(M / Mhat)))
}