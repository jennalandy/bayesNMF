#' Sample parameter: E for nth signature
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param from_prior boolean, whether to sample directly from prior distributions
#' 
#' @return vector, sampled value of E for nth signature
#' @noRd
sample_En <- function(self, n, from_prior = FALSE) {
  if (from_prior | self$params$A[1, n] == 0) {
    if (self$specs$prior == "truncnormal") {
      sampled <- truncnorm::rtruncnorm(
        1,
        mean = self$prior_params$Mu_e[n, ],
        sd = sqrt(self$prior_params$Sigmasq_e[n, ]),
        a = 0, b = Inf
      )
    } else if (self$specs$prior == "exponential") {
      sampled <- stats::rexp(self$dims$G, self$prior_params$Lambda_e[n, ])
    } else if (self$specs$prior == "gamma") {
      sampled <- stats::rgamma(
        self$dims$G,
        self$prior_params$Alpha_e[n, ], 
        self$prior_params$Beta_e[n, ]
      )
    }
    return(sampled)
  }

  if (self$specs$likelihood == "normal") {
    sampled <- sample_En_normal(self, n)
  } else if (self$specs$likelihood == "poisson" & self$specs$MH) {
    self$log("in sample_En, likelihood = poisson & MH", verbosity = 4)
    proposal <- sample_En_normal(self, n, as_proposal = TRUE)
    sampled <- MH_En_poisson(self, proposal, n)
  } else if (self$specs$likelihood == "poisson") {
    sampled <- sample_En_poisson(self, n)
  }
  return(sampled)
}

#' Sample parameter: E for nth signature (normal)
#' @description
#' Helper function to sample_En
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param as_proposal boolean, whether to sample as a proposal for Metropolis-Hastings. This affects the variance plugged in for sigmasq_kg.
#' 
#' @return vector, sampled value of E for nth signature
#' @noRd
sample_En_normal <- function(self, n, as_proposal = FALSE) {
  # Sample directly from prior when Ann = 0 or gamma = 0
  if (self$params$A[1, n] == 0 || all(self$params$P[, n] == 0)) {
    self$log("in sample_En_normal, sampling from prior", verbosity = 4)
    if (self$specs$prior == "truncnormal") {
      return(truncnorm::rtruncnorm(
        self$dims$G,
        mean = self$prior_params$Mu_e[n, ],
        sd   = sqrt(self$prior_params$Sigmasq_e[n, ]),
        a = 0, b = Inf
      ))
    } else if (self$specs$prior == "exponential") {
      return(stats::rexp(self$dims$G, self$prior_params$Lambda_e[n, ]))
    } else if (self$specs$prior == "gamma") {
      return(stats::rgamma(
        self$dims$G,
        shape = self$prior_params$Alpha_e[n, ],
        rate  = self$prior_params$Beta_e[n, ]
      ))
    }
  }

  # Otherwise sample from full conditional
  mu_sigmasq_E <- get_mu_sigmasq_En_normal(self, n, as_proposal)
  sampled <- truncnorm::rtruncnorm(
    1,
    mean = mu_sigmasq_E$mu,
    sd = sqrt(mu_sigmasq_E$sigmasq),
    a = 0,
    b = Inf
  )
  return(sampled)
}

#' Sample parameter: E for nth signature (Poisson)
#' @description
#' Helper function to sample_En
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' 
#' @return vector, sampled value of E for nth signature
#' @noRd
sample_En_poisson <- function(self, n) {
  if (self$specs$prior == "gamma") {
    shape <- sapply(1:self$dims$G, function(g) {
      self$prior_params$Alpha_e[n,g] + sum(self$params$Z[, n, g])
    })
    rate <- sapply(1:self$dims$G, function(g) {
      self$prior_params$Beta_e[n,g] +
        self$params$A[1, n] * sum(self$params$P[, n])
    })
  } else if (self$specs$prior == "exponential") {
    shape <- sapply(1:self$dims$G, function(g) {
      1 + sum(self$params$Z[, n, g])
    })
    rate <- sapply(1:self$dims$G, function(g) {
      self$prior_params$Lambda_e[n, g] +
        self$params$A[1, n] * sum(self$params$P[, n])
    })
  }
  sampled <- sapply(1:self$dims$G, function(g) {
    rgamma(1, shape[g], rate[g])
  })
  return(sampled)
}

#' Compute mean and variance for E for nth signature (normal)
#' @description
#' Helper function to sample_En_normal
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param as_proposal boolean, whether to sample as a proposal for Metropolis-Hastings. This affects the variance plugged in for sigmasq_kg.
#' 
#' @return list with mean and variance of E for nth signature
#' @noRd
get_mu_sigmasq_En_normal <- function(self, n, as_proposal = FALSE) {
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

  # broadcast Pn to residuals (M - Mhat) / sigmasq
  mu_num_term_1 <- broadcast(
    self$params$P[, n],
    across = 'columns',
    of = (self$data - Mhat_no_n) / sigmasq,
    with = '*'
  )  %>% # dim KxG
    colSums() # length G

  denom <- broadcast(
    self$params$A[1,n] * self$params$P[, n] ** 2,
    across = 'columns',
    of = 1/sigmasq,
    with = '*'
  ) %>% # dim KxG
    colSums() # length G

  if (self$specs$prior == 'exponential') {
    mu_num_term_2 <- self$prior_params$Lambda_e[n, ] # length G
    mu_E <- (mu_num_term_1 - mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom # length G
  } else if (self$specs$prior == 'truncnormal') {
    mu_num_term_2 <- self$prior_params$Mu_e[n, ] /
                     self$prior_params$Sigmasq_e[n, ] # length G
    denom <- denom + (1/self$prior_params$Sigmasq_e[n, ]) # length G
    mu_E <- (mu_num_term_1 + mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom # length G
  }

  return(list(
    mu = mu_E,
    sigmasq = sigmasq_E
  ))
}

#' Metropolis-Hastings update for E for nth signature (Poisson)
#' @description
#' Helper function to sample_En
#' 
#' @param self bayesNMF_sampler object
#' @param proposal_En vector, proposed value of E for nth signature
#' @param n signature index
#' 
#' @return vector, sampled value of E for nth signature
#' @noRd
MH_En_poisson <- function(self, proposal_En, n) {
  # For warmup samples, accept all
  if (!self$state$converged) {
    self$acceptance_rates$E_acceptance_rate[n, ] <- 1
    return(proposal_En)
  }

  proposal_E_copy <- self$params$E
  proposal_E_copy[n, ] <- proposal_En

  Mhat <- self$get_Mhat()
  proposal_Mhat <- self$get_Mhat(E = proposal_E_copy)

  # compute log likelihoods under proposal and target distributions
  loglik_poisson_old <- self$get_loglik(
    return_matrix = TRUE
  ) %>% colSums() # G-dimensional vector
  loglik_poisson_new <- self$get_loglik(
    E = proposal_E_copy, return_matrix = TRUE
  ) %>% colSums() # G-dimensional vector
  loglik_normal_old <- self$get_loglik(
    sigmasq = pmax(proposal_Mhat, 1), likelihood = "normal", return_matrix = TRUE
  ) %>% colSums() # G-dimensional vector
  loglik_normal_new <- self$get_loglik(
    E = proposal_E_copy,
    sigmasq = pmax(Mhat, 1),
    likelihood = "normal",
    return_matrix = TRUE
  ) %>% colSums() # G-dimensional vector

  # compute acceptance ratio
  # priors cancel out in ratio, same in proposal and target
  accept_ratio <- exp(
    loglik_poisson_new + loglik_normal_old -
    (loglik_poisson_old + loglik_normal_new)
  ) # G-dimensional vector
  accept_ratio <- pmin(accept_ratio, 1)
  self$acceptance_rates$E_acceptance_rate[n, ] <- accept_ratio

  # accept or reject each element of En
  u <- runif(self$dims$G)
  sampled <- self$params$E[n, ]
  sampled[u < accept_ratio] <- proposal_En[u < accept_ratio]

  return(sampled)
}