#' Sample parameter: P for nth signature
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param from_prior boolean, whether to sample directly from prior distributions
#' 
#' @return vector, sampled value of P for nth signature
#' @noRd
sample_Pn <- function(self, n, from_prior = FALSE) {
  if (from_prior | self$params$A[1, n] == 0) {
    if (self$specs$prior == "truncnormal") {
      sampled <- truncnorm::rtruncnorm(
        1,
        mean = self$prior_params$Mu_p[, n],
        sd = sqrt(self$prior_params$Sigmasq_p[, n]),
        a = 0, b = Inf
      )
    } else if (self$specs$prior == "exponential") {
      sampled <- stats::rexp(self$dims$K, self$prior_params$Lambda_p[, n])
    } else if (self$specs$prior == "gamma") {
      sampled <- stats::rgamma(
        self$dims$K,
        self$prior_params$Alpha_p[, n],
        self$prior_params$Beta_p[, n]
      )
    }
    return(sampled)
  }

  if (self$specs$likelihood == "normal") {
    sampled <- sample_Pn_normal(self, n)
  } else if (self$specs$likelihood == "poisson" & self$specs$MH) {
    self$log("in sample_Pn, likelihood = poisson & MH", verbosity = 4)
    proposal <- sample_Pn_normal(self, n, as_proposal = TRUE)
    sampled <- MH_Pn_poisson(self, proposal, n)
  } else if (self$specs$likelihood == "poisson") {
    sampled <- sample_Pn_poisson(self, n)
  }
  return(sampled)
}

#' Sample parameter: P for nth signature (normal)
#' @description
#' Helper function to sample_Pn
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param as_proposal boolean, whether to sample as a proposal for Metropolis-Hastings. This affects the variance plugged in for sigmasq_kg.
#' 
#' @return vector, sampled value of P for nth signature
#' @noRd
sample_Pn_normal <- function(self, n, as_proposal = FALSE) {
  # Sample directly from prior when Ann = 0 or gamma = 0
  if (self$params$A[1, n] == 0 || all(self$params$E[n, ] == 0)) {
    self$log("in sample_Pn_normal, sampling from prior", verbosity = 4)
    if (self$specs$prior == "truncnormal") {
      return(truncnorm::rtruncnorm(
        self$dims$K,
        mean = self$prior_params$Mu_p[, n],
        sd   = sqrt(self$prior_params$Sigmasq_p[, n]),
        a = 0, b = Inf
      ))
    } else if (self$specs$prior == "exponential") {
      return(stats::rexp(self$dims$K, self$prior_params$Lambda_p[, n]))
    } else if (self$specs$prior == "gamma") {
      return(stats::rgamma(
        self$dims$K,
        shape = self$prior_params$Alpha_p[, n],
        rate  = self$prior_params$Beta_p[, n]
      ))
    }
  }

  # Otherwise sample from full conditional
  self$log("in sample_Pn_normal, sampling from full conditional", verbosity = 4)
  mu_sigmasq_P <- get_mu_sigmasq_Pn_normal(self, n, as_proposal)
  sampled <- truncnorm::rtruncnorm(
    1,
    mean = mu_sigmasq_P$mu,
    sd = sqrt(mu_sigmasq_P$sigmasq),
    a = 0,
    b = Inf
  )
  return(sampled)
}

#' Sample parameter: P for nth signature (Poisson)
#' @description
#' Helper function to sample_Pn
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' 
#' @return vector, sampled value of P for nth signature
#' @noRd
sample_Pn_poisson <- function(self, n) {
  if (self$specs$prior == "gamma") {
    shape <- sapply(1:self$dims$K, function(k) {
      self$prior_params$Alpha_p[k, n] + sum(self$params$Z[k, n, ])
    })
    rate <- sapply(1:self$dims$K, function(k) {
      self$prior_params$Beta_p[k, n] +
        self$params$A[1, n] * sum(self$params$E[n, ])
    })
  } else if (self$specs$prior == "exponential") {
    shape <- sapply(1:self$dims$K, function(k) {
      1 + sum(self$params$Z[k,n,])
    })
    rate <- sapply(1:self$dims$K, function(k) {
      self$prior_params$Lambda_p[k, n] +
        self$params$A[1, n] * sum(self$params$E[n, ])
    })
  }
  sampled <- sapply(1:self$dims$K, function(k) {
    rgamma(1, shape[k], rate[k])
  })
  return(sampled)
}

#' Compute mean and variance for P for nth signature (normal)
#' @description
#' Helper function to sample_Pn_normal
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param as_proposal boolean, whether to sample as a proposal for Metropolis-Hastings. This affects the variance plugged in for sigmasq_kg.
#' 
#' @return list with mean and variance of P for nth signature
#' @noRd
get_mu_sigmasq_Pn_normal <- function(self, n, as_proposal = FALSE) {
  self$log("in get_mu_sigmasq_Pn_normal", verbosity = 4)

  # compute Mhat
  Mhat <- self$get_Mhat()
  if (as_proposal) {
    # in proposal, mean = variance to match Poisson assumption
    sigmasq <- Mhat
  } else {
    # otherwise, create sigmasq matrix with each row equal to self$params$sigmasq
    sigmasq <- matrix(
      rep(self$params$sigmasq, times = self$dims$K),
      nrow = self$dims$K,
      byrow = TRUE
    )
  }

  # compute Mhat excluding nth signature
  # make a copy of the current A matix and put 0 in the nth column
  A_copy <- self$params$A; A_copy[1, n] <- 0
  Mhat_no_n <- self$get_Mhat(A = A_copy)

  # broadcast En to residuals (M - Mhat) / sigmasq
  mu_num_term_1 <- broadcast(
    self$params$E[n, ],
    across = 'rows',
    of = (self$data - Mhat_no_n) / sigmasq,
    with = '*'
  )  %>% # dim KxG
    rowSums() # length K

  denom <- broadcast(
    self$params$A[1,n] * self$params$E[n, ] ** 2,
    across = 'rows',
    of = 1/sigmasq,
    with = '*'
  ) %>% # dim KxG
    rowSums() # length K

  if (self$specs$prior == 'exponential') {
    mu_num_term_2 <- self$prior_params$Lambda_p[, n] # length K
    mu_P <- (mu_num_term_1 - mu_num_term_2) / denom # length K
    sigmasq_P <- 1 / denom # length K
  } else if (self$specs$prior == 'truncnormal') {
    mu_num_term_2 <- self$prior_params$Mu_p[, n] / 
                     self$prior_params$Sigmasq_p[,n] # length K
    denom <- denom + (1 / self$prior_params$Sigmasq_p[,n]) # length K
    mu_P <- (mu_num_term_1 + mu_num_term_2) / denom # length K
    sigmasq_P <- 1 / denom # length K
  }

  return(list(
    mu = mu_P,
    sigmasq = sigmasq_P
  ))
}

#' Metropolis-Hastings update for P for nth signature (Poisson)
#' @description
#' Helper function to sample_Pn
#' 
#' @param self bayesNMF_sampler object
#' @param proposal_Pn vector, proposed value of P for nth signature
#' @param n signature index
#' 
#' @return vector, sampled value of P for nth signature
#' @noRd
MH_Pn_poisson <- function(self, proposal_Pn, n) {
  # for warmup samples, accept all
  if (!self$state$converged) {
    self$acceptance_rates$P_acceptance_rate[, n] <- 1
    return(proposal_Pn)
  }

  proposal_P_copy <- self$params$P
  proposal_P_copy[, n] <- proposal_Pn

  Mhat <- self$get_Mhat()
  proposal_Mhat <- self$get_Mhat(P = proposal_P_copy)

  # compute log likelihoods under proposal and target distributions
  loglik_poisson_old <- self$get_loglik(
    return_matrix = TRUE
  ) %>% rowSums() # K-dimensional vector
  loglik_poisson_new <- self$get_loglik(
    P = proposal_P_copy, return_matrix = TRUE
  ) %>% rowSums() # K-dimensional vector
  loglik_normal_old <- self$get_loglik(
    # conditional on proposed P, which goes into sigmasq
    sigmasq = pmax(proposal_Mhat, 1),
    likelihood = "normal", 
    return_matrix = TRUE
  ) %>% rowSums() # K-dimensional vector
  loglik_normal_new <- self$get_loglik(
    P = proposal_P_copy,
    # conditional on previous P, which goes into sigmasq
    sigmasq = pmax(Mhat, 1),
    likelihood = "normal",
    return_matrix = TRUE
  ) %>% rowSums() # K-dimensional vector

  # compute acceptance ratio
  # priors cancel out in ratio, same in proposal and target
  accept_ratio <- exp(
    loglik_poisson_new + loglik_normal_old -
    (loglik_poisson_old + loglik_normal_new)
  ) # K-dimensional vector
  accept_ratio <- pmin(accept_ratio, 1)
  self$acceptance_rates$P_acceptance_rate[, n] <- accept_ratio

  # accept or reject each element of Pn
  u <- runif(self$dims$K)
  sampled <- self$params$P[, n]
  sampled[u < accept_ratio] <- proposal_Pn[u < accept_ratio]

  return(sampled)
}