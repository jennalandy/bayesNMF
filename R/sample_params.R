# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PRIVATE METHODS ##############
####################################

#' Initialize parameters
#' @description
#' Initialize parameters based on specified initial values or sampled from prior distributions if not provided
#' 
#' @param self bayesNMF_sampler object
#' @param init_params list, initial values for parameters
#' 
#' @return None, updates self$params
#' @noRd
init_params_ <- function(self, init_params) {
  self$params <- init_params

  if (!("P" %in% names(self$params))) {
    self$params$P <- matrix(nrow = self$dims$K, ncol = self$dims$N)
  }
  if (!("E" %in% names(self$params))) {
    self$params$E <- matrix(nrow = self$dims$N, ncol = self$dims$G)
  }
  if (self$specs$likelihood == "poisson" & !self$specs$MH) {
    if (!("Z" %in% names(self$params))) {
      self$params$Z <- array(dim = c(self$dims$K, self$dims$N, self$dims$G))
    }
  }
  if (self$specs$likelihood == "normal") {
    if (!("sigmasq" %in% names(self$params))) {
      self$params$sigmasq <- vector(mode = "numeric", length = self$dims$G)
    }
  }
  if (!("A" %in% names(self$params))) {
    self$params$A <- matrix(nrow = 1, ncol = self$dims$N)
  }
  if (!("R" %in% names(self$params))) {
    self$params$R <- self$dims$N
  }
}

#' Sample parameters from full conditionals (one iteration of Gibbs sampling)
#' 
#' @param self bayesNMF_sampler object
#' @param skip list, names of parameters to skip
#' @param from_prior boolean, whether to sample directlyfrom prior distributions
#' 
#' @return None, updates self$params
#' @noRd
sample_params_ <- function(self, skip, from_prior = FALSE) {
  updates <- get_updates_()

  if (!('P' %in% skip)) {
    self$log("Sampling P", verbosity = 3)
    for (n in 1:self$dims$N) {
      self$params$P[, n] <- updates$Pn(self, n, from_prior)
    }
  }
  if (!('E' %in% skip)) {
    self$log("Sampling E", verbosity = 3)
    for (n in 1:self$dims$N) {
      self$params$E[n, ] <- updates$En(self, n, from_prior)
    }
  }

  if (!('A' %in% skip) & self$specs$learning_rank) {
    self$log("Sampling R", verbosity = 3)
    self$params$R <- updates$R(self, from_prior)

    self$log("Sampling A", verbosity = 3)
    for (n in 1:self$dims$N) {
      self$params$A[1, n] <- updates$An(self, n, from_prior)
    }
  } else if (!('A' %in% skip) & !self$specs$learning_rank & from_prior) {
    self$params$A <- matrix(1, nrow = 1, ncol = self$dims$N)
  }

  if (self$specs$likelihood == 'poisson' & !self$specs$MH & !('Z' %in% skip)) {
    for (k in 1:self$dims$K) {
      for (g in 1:self$dims$G) {
        self$params$Z[k, , g] <- updates$Zkg(self, k, g)
      }
    }
  }
  if (self$specs$likelihood == 'normal' & !('sigmasq' %in% skip)) {
    self$params$sigmasq <- updates$sigmasq(self)
  }
}

#' Sample parameter: A for nth signature
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' @param n signature index
#' @param from_prior boolean, whether to sample directly from prior distributions
#' 
#' @return integer, sampled value of A for nth signature
#' @noRd
sample_An <- function(self, n, from_prior = FALSE) {
  prior_prob_1 <- compute_prior_prob_1(self$params$R, self$dims$N)
  if (from_prior) {
    sampled <- rbinom(1, 1, prior_prob_1)
    return(sampled)
  }

  temperature <- self$temperature_schedule[self$state$iter]

  # make copies of A with 1 and 0 for signature n
  A0 <- self$params$A; A0[1, n] <- 0
  A1 <- self$params$A; A1[1, n] <- 1

  # compute log likelihoods
  loglik_0 <- self$get_loglik(A = A0)
  loglik_1 <- self$get_loglik(A = A1)

  if (self$specs$rank_method == "SBFI") {
    n_params_0 <- sum(A0) * (self$dims$G + self$dims$K)
    n_params_1 <- sum(A1) * (self$dims$G + self$dims$K)

    neg_half_BIC_0 <- loglik_0 - n_params_0 * log(self$dims$G) / 2
    neg_half_BIC_1 <- loglik_1 - n_params_1 * log(self$dims$G) / 2

    log_p0 = log(1 - prior_prob_1) + temperature * neg_half_BIC_0
    log_p1 = log(prior_prob_1) + temperature * neg_half_BIC_1
  } else {
    log_p0 = log(1 - prior_prob_1) + temperature * loglik_0
    log_p1 = log(prior_prob_1) + temperature * loglik_1
  }

  log_p = log_p1 - sumLog(c(log_p0, log_p1))
  p = exp(log_p)

  # handle overflow
  if (is.na(p)) {
    self$log("Overflow in sample_An", verbosity = 1)
    self$log(paste0("Previous A: ", paste(self$params$A, collapse = ', ')), verbosity = 1)
    self$log(paste0("n: ", n), verbosity = 1)
    self$log(paste0("log_p1: ", log_p1), verbosity = 1)
    self$log(paste0("log_p0: ", log_p0), verbosity = 1)
    self$log(paste0("log_p: ", log_p), verbosity = 1)
    self$log(paste0("p: ", p), verbosity = 1)
    self$log(paste0("prior_prob_1: ", prior_prob_1), verbosity = 1)
    self$log(paste0("temperature: ", temperature), verbosity = 1)
    self$log(paste0("loglik_0: ", loglik_0), verbosity = 1)
    self$log(paste0("loglik_1: ", loglik_1), verbosity = 1)
    if (is.na(log_p1) & is.na(log_p0)) {
      p = 0.5
    } else if (is.na(log_p1)) {
      p = 0
    } else if (is.na(log_p0)) {
      p = 1
    } else if (log_p1 > log_p0) {
      p = 1
    } else if (log_p1 < logp0) {
      p = 0
    } else {
      p = 0.5
    }
    self$log(paste0("Resolved p is ", p), verbosity = 1)
  }

  sampled <- rbinom(1, 1, p)
  return(sampled)
}

#' Compute prior probability of inclusion
#' @description
#' Helper function to sample_An
#' 
#' @param R integer, expected rank
#' @param N integer, number of samples
#' @param clip_val double, clipping value (default 0.4)
#' 
#' @return double, prior probability of inclusion
#' @noRd
compute_prior_prob_1 <- function(R, N, clip_val = 0.4) {
  # prior probability of inclusion is the same for each signature
  # defined by R, which is expected rank
  prior_prob_1 <- R / N

  # clip prior_prob_1 away from 0 and 1 by clip_val/N
  prior_prob_1 <- max(prior_prob_1, clip_val / N)
  prior_prob_1 <- min(prior_prob_1, 1 - clip_val / N)
  return(prior_prob_1)
}

#' Compute log sum of exponential
#' @description
#' Given a vector of log values log(x_i), i = 1, ..., n, 
#' return the log of the sum of values log(x_1 + x_2 + ... + x_n)
#' Helper function to sample_An
#' 
#' @param vec vector, values to compute log sum of exponential
#' 
#' @return double, log sum of exponential
#' @noRd
sumLog <- function(vec) {
  ord <- sort(vec, decreasing = TRUE)
  s <- ord[1]
  for (i in 2:length(vec)) {
    s <- s + log( 1 + exp(ord[i] - s))
  }
  return(s)
}

#' Sample parameter: R (expected rank)
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' @param from_prior boolean, whether to sample directly from prior distributions
#' 
#' @return integer, sampled value of R
#' @noRd
sample_R <- function(self, from_prior = FALSE) {
  range_R = 0:self$dims$N
  temperature <- self$temperature_schedule[self$state$iter]

  # compute posterior probability of each expected rank
  probs <- sapply(range_R, function(r) {
    prior_prob_R <- 1/length(range_R)
    if (from_prior) {
      return(prior_prob_R)
    }

    prob_A1 <- compute_prior_prob_1(r, self$dims$N)
    likelihood_A <- (
      prob_A1 ** sum(self$params$A) *
        (1 - prob_A1) ** (self$dims$N - sum(self$params$A))
    ) ** temperature
    posterior_prob_R <- prior_prob_R * likelihood_A
    return(posterior_prob_R)
  })

  # normalize probabilities and sample new expected rank
  probs <- probs / sum(probs)
  sampled <- sample(range_R, size = 1, prob = probs)
  return(sampled)
}

#' Sample parameter: Z for kth signature and gth genome
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' @param k signature index
#' @param g genome index
#' 
#' @return integer, sampled value of Z for kth signature and gth genome
#' @noRd
sample_Zkg <- function(self, k, g) {
  probs <- sapply(1:self$dims$N, function(n) {
    self$params$P[k, n] * self$params$A[1, n] * self$params$E[n, g]
  })
  if (sum(probs) == 0) {
    # happens if all factors are excluded
    # setting Zs to 0 means P, E sampled from their priors
    return(rep(0, self$dims$N))
  }
  probs <- probs / sum(probs)
  sampled <- rmultinom(1, size = self$data[k, g], prob = probs)
  return(sampled)
}

#' Sample parameter: sigmasq for gth genome
#' @description
#' Helper function to sample_params_
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return vector, sampled value of sigmasq for gth genome
#' @noRd
sample_sigmasq <- function(self) {
  Mhat <- self$get_Mhat()
  sigmasq <- sapply(1:self$dims$G, function(g) {
    residuals <- self$data[, g] - Mhat[, g]
    invgamma::rinvgamma(
      1,
      shape = self$prior_params$Alpha[g] + self$dims$K / 2,
      rate = self$prior_params$Beta[g] + (1 / 2) * sum(residuals ** 2)
    )
  })
  return(sigmasq)
}

#' Get update functions for all parameters
#' @description
#' Helper function to sample_params_
#' 
#' @return list, update functions for all parameters
#' @noRd
get_updates_ <- function() {
  return(list(
    Pn = sample_Pn,
    En = sample_En,
    An = sample_An,
    R = sample_R,
    Zkg = sample_Zkg,
    sigmasq = sample_sigmasq
  ))
}