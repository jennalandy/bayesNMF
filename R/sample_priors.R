# These functions are copied over to be methods of the bayesNMF_sampler class

####################################
###### PRIVATE METHODS ##############
####################################

#' Initialize prior parameters
#' @description
#' Initialize prior parameters based on specified initial values or sampled from hyperprior distributions if not provided
#' 
#' @param self bayesNMF_sampler object
#' 
#' @return None, updates self$prior_params
#' @noRd
init_prior_params_ <- function(self, init_prior_params) {
  self$prior_params <- init_prior_params

  if (self$specs$prior == "truncnormal") {
    if (!("Mu_p" %in% names(self$prior_params))) {
      self$prior_params$Mu_p <- matrix(nrow = self$dims$K, ncol = self$dims$N)
    }
    if (!("Sigmasq_p" %in% names(self$prior_params))) {
      self$prior_params$Sigmasq_p <- matrix(nrow = self$dims$K, ncol = self$dims$N)
    }
    if (!("Mu_e" %in% names(self$prior_params))) {
      self$prior_params$Mu_e <- matrix(nrow = self$dims$N, ncol = self$dims$G)
    }
    if (!("Sigmasq_e" %in% names(self$prior_params))) {
      self$prior_params$Sigmasq_e <- matrix(nrow = self$dims$N, ncol = self$dims$G)
    }

    for (n in 1:self$dims$N) {
      if (any(is.na(self$prior_params$Mu_p[, n]))) {
        self$prior_params$Mu_p[, n] <- rnorm(
          self$dims$K,
          self$hyperprior_params$M_p[, n],
          sqrt(self$hyperprior_params$S_p[, n])
        )
      }
      if (any(is.na(self$prior_params$Sigmasq_p[, n]))) {
        self$prior_params$Sigmasq_p[, n] <- invgamma::rinvgamma(
          n = self$dims$K,
          shape = self$hyperprior_params$A_p[, n],
          rate = self$hyperprior_params$B_p[, n]
        )
      }
      if (any(is.na(self$prior_params$Mu_e[n, ]))) {
        self$prior_params$Mu_e[n, ] <- rnorm(
          self$dims$G,
          self$hyperprior_params$M_e[n, ],
          sqrt(self$hyperprior_params$S_e[n, ])
        )
      }
      if (any(is.na(self$prior_params$Sigmasq_e[n, ]))) {
        self$prior_params$Sigmasq_e[n, ] <- invgamma::rinvgamma(
          n = self$dims$G,
          shape = self$hyperprior_params$A_e[n, ],
          rate = self$hyperprior_params$B_e[n, ]
        )
      }
    }
  } else if (self$specs$prior == "exponential") {
    if (!("Lambda_p" %in% names(self$prior_params))) {
      self$prior_params$Lambda_p <- matrix(nrow = self$dims$K, ncol = self$dims$N)
    }
    if (!("Lambda_e" %in% names(self$prior_params))) {
      self$prior_params$Lambda_e <- matrix(nrow = self$dims$N, ncol = self$dims$G)
    }
    
    for (n in 1:self$dims$N) {
      if (any(is.na(self$prior_params$Lambda_p[, n]))) {
        self$prior_params$Lambda_p[, n] <- rgamma(
          n = self$dims$K,
          shape = self$hyperprior_params$A_p[, n],
          rate = self$hyperprior_params$B_p[, n]
        )
      }
      if (any(is.na(self$prior_params$Lambda_e[n, ]))) {
        self$prior_params$Lambda_e[n, ] <- rgamma(
          n = self$dims$G,
          shape = self$hyperprior_params$A_e[n, ],
          rate = self$hyperprior_params$B_e[n, ]
        )
      }
    }
  } else if (self$specs$prior == "gamma") {
    if (!("Alpha_p" %in% names(self$prior_params))) {
      self$prior_params$Alpha_p <- matrix(nrow = self$dims$K, ncol = self$dims$N)
    }
    if (!("Beta_p" %in% names(self$prior_params))) {
      self$prior_params$Beta_p <- matrix(nrow = self$dims$K, ncol = self$dims$N)
    }
    
    if (!("Alpha_e" %in% names(self$prior_params))) {
      self$prior_params$Alpha_e <- matrix(nrow = self$dims$N, ncol = self$dims$G)
    }
    if (!("Beta_e" %in% names(self$prior_params))) {
      self$prior_params$Beta_e <- matrix(nrow = self$dims$N, ncol = self$dims$G)
    }

    for (n in 1:self$dims$N) {
      if (any(is.na(self$prior_params$Beta_p[, n]))) {
        self$prior_params$Beta_p[, n] <- rgamma(
          n = self$dims$K,
          shape = self$hyperprior_params$A_p[, n],
          rate = self$hyperprior_params$B_p[, n]
        )
      }
      if (any(is.na(self$prior_params$Alpha_p[, n]))) {
        self$prior_params$Alpha_p[, n] <- rgamma(
          n = self$dims$K,
          shape = self$hyperprior_params$C_p[, n],
          rate = self$hyperprior_params$D_p[, n]
        )
      }
      if (any(is.na(self$prior_params$Beta_e[n, ]))) {
        self$prior_params$Beta_e[n, ] <- rgamma(
          n = self$dims$G,
          shape = self$hyperprior_params$A_e[n, ],
          rate = self$hyperprior_params$B_e[n, ]
        )
      }
      if (any(is.na(self$prior_params$Alpha_e[n, ]))) {
        self$prior_params$Alpha_e[n, ] <- rgamma(
          n = self$dims$G,
          shape = self$hyperprior_params$C_e[n, ],
          rate = self$hyperprior_params$D_e[n, ]
        )
      }
    }
  }

  if (self$specs$likelihood == "normal") {
    if (!("Alpha" %in% self$prior_params)) {
      self$prior_params$Alpha <- rep(self$prior_params$alpha, length = self$dims$G)
    }
    if (!("Beta" %in% self$prior_params)) {
      self$prior_params$Beta <- rep(self$prior_params$beta, length = self$dims$G)
    }
  }
}

#' Sample prior parameters from full conditionals
#' 
#' @param self bayesNMF_sampler object
#' @param skip list, names of prior parameters to skip
#' 
#' @return None, updates self$prior_params
#' @noRd
sample_prior_params_ <- function(self, skip) {
  if (self$specs$prior == "truncnormal") {
    for (n in 1:self$dims$N) {
      if (!("Mu_p" %in% skip)) {
        self$log("sampling Mu_p", verbosity = 3)
        self$prior_params$Mu_p[, n] <- sample_Mu_Pn(self, n)
      }
      if (!("Mu_e" %in% skip)) {
        self$log("sampling Mu_e", verbosity = 3)
        self$prior_params$Mu_e[n, ] <- sample_Mu_En(self, n)
      }
      if (!("Sigmasq_p" %in% skip)) {
        self$log("sampling Sigmasq_p", verbosity = 3)
        self$prior_params$Sigmasq_p[, n] <- sample_Sigmasq_Pn(self, n)
      }
      if (!("Sigmasq_e" %in% skip)) {
        self$log("sampling Sigmasq_e", verbosity = 3)
        self$prior_params$Sigmasq_e[n, ] <- sample_Sigmasq_En(self, n)
      }
    }
  } else if (self$specs$prior == "exponential") {
    for (n in 1:self$dims$N) {
      if (!("Lambda_p" %in% skip)) {
        self$prior_params$Lambda_p[, n] <- sample_Lambda_Pn(self, n)
      }
      if (!("Lambda_e" %in% skip)) {
        self$prior_params$Lambda_e[n, ] <- sample_Lambda_En(self, n)
      }
    }
  } else if (self$specs$prior == "gamma") {
    for (n in 1:self$dims$N) {
      if (!("Beta_p" %in% skip)) {
        self$prior_params$Beta_p[, n] <- sample_Beta_Pn(self, n)
      }
      for (k in 1:self$dims$K) {
        if (!("Alpha_p" %in% skip)) {
          self$prior_params$Alpha_p[k, n] <- sample_Alpha_Pkn(self, k, n)
        }
      }

      if (!("Beta_e" %in% skip)) {
        self$prior_params$Beta_e[n, ] <- sample_Beta_En(self, n)
      }
      for (g in 1:self$dims$G) {
        if (!("Alpha_e" %in% skip)) {
          self$prior_params$Alpha_e[n,g] <- sample_Alpha_Eng(self, n, g)
        }
      }
    }
  }
}

#####################################################################
####################### Truncated Normal Prior ######################
#####################################################################

#' Sample Prior Parameter: Mu for Pn
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length K
#' @noRd
sample_Mu_Pn <- function(self, n) {
    num <- self$hyperprior_params$M_p[, n]/self$hyperprior_params$S_p[, n] + 
           self$params$P[, n]/self$prior_params$Sigmasq_p[, n]
    denom <- 1/self$hyperprior_params$S_p[, n] + 
             1/self$prior_params$Sigmasq_p[, n]
    rnorm(self$dims$K, num/denom, 1/denom)
}

#' Sample Prior Parameter: Mu for En
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length G
#' @noRd
sample_Mu_En <- function(self, n) {
    num <- self$hyperprior_params$M_e[n, ]/self$hyperprior_params$S_e[n, ] +
           self$params$E[n, ]/self$prior_params$Sigmasq_e[n, ]
    denom <- 1/self$hyperprior_params$S_e[n, ] +
            1/self$prior_params$Sigmasq_e[n, ]
    rnorm(self$dims$G, num/denom, 1/denom)
}

#' Sample Prior Parameter: Sigmasq for Pn
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length K
#' @noRd
sample_Sigmasq_Pn <- function(self, n) {
    invgamma::rinvgamma(
        n = self$dims$K,
        shape = self$hyperprior_params$A_p[, n] + 1/2,
        rate = self$hyperprior_params$B_p[, n] +
               (self$params$P[, n] - self$prior_params$Mu_p[, n])**2/2
    )
}

#' Sample Prior Parameter: Sigmasq for En
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length G
#' @noRd
sample_Sigmasq_En <- function(self, n) {
    invgamma::rinvgamma(
        n = self$dims$G,
        shape = self$hyperprior_params$A_e[n, ] + 1/2,
        rate = self$hyperprior_params$A_e[n, ] +
               (self$params$E[n, ] - self$prior_params$Mu_e[n, ])**2/2
    )
}

#####################################################################
####################### Exponential Prior ######################
#####################################################################

#' Sample Prior Parameter: Lambda for Pn
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length K
#' @noRd
sample_Lambda_Pn <- function(self, n) {
  sampled <- rgamma(
    self$dims$K,
    self$hyperprior_params$A_p[, n] + 1,
    self$hyperprior_params$B_p[, n] + self$params$P[, n]
  )
  return(sampled)
}

#' Sample Prior Parameter: Lambda for En
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length G
#' @noRd
sample_Lambda_En <- function(self, n) {
    sampled <- rgamma(
      self$dims$G,
      self$hyperprior_params$A_e[n, ] + 1,
      self$hyperprior_params$B_e[n, ] + self$params$E[n, ]
    )
    return(sampled)
}


#####################################################################
####################### Gamma Prior ######################
#####################################################################

#' Sample Prior Parameter: Beta for Pn
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length K
#' @noRd
sample_Beta_Pn <- function(self, n) {
  rgamma(
    self$dims$K,
    self$hyperprior_params$A_p[, n] + self$prior_params$Alpha_p[, n],
    self$hyperprior_params$B_p[, n] + self$params$P[, n]
  )
}

#' Sample Prior Parameter: Beta for En
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#'
#' @return vector of length G
#' @noRd
sample_Beta_En <- function(self, n) {
  rgamma(
    self$dims$G,
    self$hyperprior_params$A_e[n, ] + self$prior_params$Alpha_e[n, ],
    self$hyperprior_params$B_e[n, ] + self$params$E[n, ]
  )
}

#' Sample Prior Parameter: Alpha for Pkn
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param k mutation type index
#' @param n signature index
#'
#' @return scalar
#' @noRd
sample_Alpha_Pkn <- function(self, k, n) {
  logpdf_prop <- function(x) {
    (self$hyperprior_params$C_p[k, n] - 1) * log(x) -
      self$hyperprior_params$D_p[k, n] * x +
      x * log(self$prior_params$Beta_p[k, n]) +
      (x - 1) * log(self$params$P[k, n]) -
      lgamma(x)
  }

  armspp::arms(
    n_samples = 1,
    log_pdf = logpdf_prop,
    lower = 1e-3,
    upper = 10000
  )
}

#' Sample Prior Parameter: Alpha for Eng
#' @description
#' Helper function to sample_prior_params_ 
#'
#' @param n signature index
#' @param g tumor genome index
#'
#' @return scalar
#' @noRd
sample_Alpha_Eng <- function(self, n, g) {
  logpdf_prop <- function(x) {
    (self$hyperprior_params$C_e[n,g] - 1) * log(x) -
      self$hyperprior_params$D_e[n,g] * x +
      x * log(self$prior_params$Beta_e[n,g]) +
      (x - 1) * log(self$params$E[n,g]) -
      lgamma(x)
  }

  armspp::arms(
    n_samples = 1,
    log_pdf = logpdf_prop,
    lower = 1e-3,
    upper = 10000
  )
}
