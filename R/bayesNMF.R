#' get mu and sigmasq for P[,n] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_exponential <- function(n, M, Theta, gamma = 1) {
    Mhat_no_n <- Theta$P[, -n] %*% Theta$E[-n, ]
    sum_E_sq <- gamma * sum(Theta$E[n, ] ** 2)

    # compute mean
    mu_num_term_1 <- gamma * sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$E[n, ], # length G
        "*"
    ) %>% # dim KxG
        rowSums() # length K

    mu_num_term_2 <- Theta$Lambda_p[, n] * Theta$sigmasq # length K
    mu_P <- (mu_num_term_1 - mu_num_term_2) / sum_E_sq # length K
    sigmasq_P <- Theta$sigmasq / sum_E_sq # length K

    return(list(
        mu = mu_P,
        sigmasq = sigmasq_P
    ))
}

#' get mu and sigmasq for P[,n] in Normal-TruncNormal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_truncnormal <- function(n, M, Theta, gamma = 1) {
    Mhat_no_n <- Theta$P[, -n] %*% diag(Theta$A[1, -n]) %*% Theta$E[-n, ]

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * (1/Theta$sigmasq) * (sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$E[n, ], # length G
        "*"
    ) %>% # dim KxG
        rowSums()) # length K
    mu_num_term_2 <- Theta$Mu_p[, n] / Theta$Sigmasq_p[,n] # length K
    denom <- (1/Theta$Sigmasq_p[,n]) + gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2) / Theta$sigmasq


    mu_P <- (mu_num_term_1 + mu_num_term_2) / denom # length K
    sigmasq_P <- 1 / denom # length K

    return(list(
        mu = mu_P,
        sigmasq = sigmasq_P
    ))
}

#' sample P[,n] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_normal <- function(n, M, Theta, prior = 'truncnormal', gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_truncnormal(n, M, Theta, gamma = gamma)
    } else if (prior == 'exponential') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_exponential(n, M, Theta, gamma = gamma)
    }


    mu_P = mu_sigmasq_P$mu
    sigmasq_P = mu_sigmasq_P$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_P, sd = sqrt(sigmasq_P), a = 0, b = Inf)
}

#' sample P[,n] for Poisson likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_poisson <- function(n, M, Theta, dims, prior = 'gamma', gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, Theta$Alpha_p[k,n] + gamma * sum(Theta$Z[k,n,]), Theta$Beta_p[k,n] + gamma * sum(Theta$E[n,]))
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, 1 + gamma * sum(Theta$Z[k,n,]), Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,]))
        })
    }
    return(sampled)
}

#' sample P[,n] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','exponential')
#' @param prior string, one of c('gamma','exponential','truncnormal')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn <- function(n, M, Theta, dims, likelihood = 'normal', prior = 'truncnormal', gamma = 1) {
    if (likelihood == 'normal') {
        sample_Pn_normal(n, M, Theta, prior, gamma)
    } else if (likelihood == 'poisson') {
        sample_Pn_poisson(n, M, Theta, dims, prior, gamma)
    }
}

#' Get mu and sigmasq for E[n,] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_exponential <- function(n, M, Theta, gamma = 1) {
    Mhat_no_n <- Theta$P[, -n] %*% Theta$E[-n, ]

    # compute mean
    mu_num_term_1 <- gamma * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n] / Theta$sigmasq, # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_2 <- Theta$Lambda_e[n, ] # length G
    denom <- sum(gamma * Theta$P[, n] ** 2 / Theta$sigmasq)

    mu_E <- (mu_num_term_1 - mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom

    return(list(
        mu = mu_E,
        sigmasq = sigmasq_E
    ))
}

#' Get mu and sigmasq for E[n,] in Normal-TruncNormal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_truncnormal <- function(n, M, Theta, gamma = 1) {
    Mhat_no_n <- Theta$P[, -n] %*% diag(Theta$A[1, -n]) %*% Theta$E[-n, ]

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n] / Theta$sigmasq, # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_2 <- Theta$Mu_e[n, ] / Theta$Sigmasq_e[n,] # length G
    denom <- (1/Theta$Sigmasq_e[n,]) + gamma * sum(Theta$A[1,n] * Theta$P[, n] ** 2 / Theta$sigmasq)

    mu_E <- (mu_num_term_1 + mu_num_term_2) / denom # length G
    sigmasq_E <- 1 / denom

    return(list(
        mu = mu_E,
        sigmasq = sigmasq_E
    ))
}

#' Sample E[n,] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_normal <- function(n, M, Theta, prior = 'truncnormal', gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_truncnormal(n, M, Theta, gamma = gamma)
    } else if (prior == 'exponential') {
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_exponential(n, M, Theta, gamma = gamma)
    }
    mu_E = mu_sigmasq_E$mu
    sigmasq_E = mu_sigmasq_E$sigmasq

    # sample from truncated normal
    truncnorm::rtruncnorm(1, mean = mu_E, sd = sqrt(sigmasq_E), a = 0, b = Inf)
}

#' Sample E[n,g] for Poisson likelihood
#'
#' @param n signature index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_poisson <- function(n, M, Theta, dims, prior = 'gamma', gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(1, Theta$Alpha_e[n,g] + gamma*sum(Theta$Z[,n,g]), Theta$Beta_e[n,g] + gamma*sum(Theta$P[,n]))
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(1, 1 + gamma * sum(Theta$Z[,n,g]), Theta$Lambda_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n]))
        })
    }
    return(sampled)
}


#' Sample E[n,] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En <- function(n, M, Theta, dims, likelihood = 'normal', prior = 'truncnormal', gamma = 1) {
    if (likelihood == 'normal') {
        sample_En_normal(n, M, Theta, prior = prior, gamma = gamma)
    } else if (likelihood == 'poisson') {
        sample_En_poisson(n, M, Theta, dims, prior = prior, gamma = gamma)
    }
}

#' sample sigmasq for Normal likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_sigmasq_normal <- function(M, Theta, dims, gamma = 1){
    Mhat <- Theta$P %*% diag(Theta$A[1, ]) %*% Theta$E
    sigmasq <- sapply(1:dims$K, function(k) {
        sigmasq_k <- invgamma::rinvgamma(
            n = 1,
            shape = Theta$Alpha[k] + gamma * dims$G / 2,
            scale = Theta$Beta[k] + gamma * sum(((M - Mhat)[k,])**2) / 2
        )
    })

    return(sigmasq)
}


#' sample Z[k,,g] for Poisson likelihood
#'
#' @param k mutation type index
#' @param g tumor genome index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return vector length K
#' @noRd
sample_Zkg_poisson <- function(k, g, M, Theta, dims){
    probs = sapply(1:dims$N, function(n) {
        Theta$P[k,n] * Theta$E[n,g]
    })
    probs = probs/sum(probs)
    rmultinom(1, size = M[k,g], prob = probs)
}

#' Sample An
#'
#' @param n integer, signature index
#' @param M list of matrices
#' @param Theta list of parameters
#' @param dims list of dimension values
#' @param logfac vector, logfac[i] = log(i!), only needed for `likelihood = 'poisson'`
#' @param likelihood string, one of c('normal','poisson')
#' @param gamma double, tempering parameter
#'
#' @return integer
#' @noRd
sample_An <- function(n, M, Theta, dims, logfac, likelihood = 'normal', gamma = 1) {
    Theta_A0 <- Theta
    Theta_A0$A[1,n] <- 0
    Theta_A1 <- Theta
    Theta_A1$A[1,n] <- 1

    if (likelihood == 'normal') {
        loglik_0 <- get_loglik_normal(M, Theta_A0, dims)
        loglik_1 <- get_loglik_normal(M, Theta_A1, dims)
    } else if (likelihood == 'poisson') {
        loglik_0 <- get_loglik_poisson(M, Theta_A0, dims, logfac)
        loglik_1 <- get_loglik_poisson(M, Theta_A1, dims, logfac)
    }

    log_p0 = log(1 - Theta$q[1,n]) + gamma * loglik_0
    log_p1 = log(Theta$q[1,n]) + gamma * loglik_1

    log_p = log_p1 - sumLog(c(log_p0, log_p1))
    p = exp(log_p)
    return(sample(c(0, 1), size = 1, prob = c(1-p, p)))
}

#' Sample qn
#'
#' @param n integer, signature index
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return double
#' @noRd
sample_qn <- function(n, Theta, gamma = 1) {
    rbeta(1, Theta$a + gamma*Theta$A[1, n], Theta$b - gamma*Theta$A[1, n] + 1)
}


#' get Normal log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return scalar
#' @noRd
get_loglik_normal <- function(M, Theta, dims) {
    - dims$G * sum(log(2 * pi * Theta$sigmasq)) / 2 -
        sum(sweep(
            (M - Theta$P %*% diag(Theta$A[1, ]) %*% Theta$E)**2,
            1,
            1/(2 * Theta$sigmasq), # Length K
            '*'
        ))
}


#' get Poisson log likelihood
#'
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param logfac vector, logfac[i] = log(i!)
#'
#' @return scalar
#' @noRd
get_loglik_poisson <- function(M, Theta, dims, logfac) {
    - sum((Theta$P %*% Theta$E)) +
        sum(M * log(Theta$P %*% Theta$E)) -
        sum(logfac[M])
}


#' Set prior parameters for Truncated Normal prior
#'
#' @param dims named list of dimensions N, K, G
#' @param mu_p see `Mu_p`
#' @param Mu_p mean for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `mu_p`
#' @param sigmasq_p see `Sigmasq_p`
#' @param Sigmasq_p variance for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `sigmasq_p`
#' @param mu_e see `Mu_e`
#' @param Mu_e mean for the truncated normal prior on `E`, matrix
#' size N x G. Defaults to all same value `mu_e`
#' @param sigmasq_e see `Sigmasq_e`
#' @param Sigmasq_e variance for the truncated normal prior on `P`, matrix
#' size N x G. Defaults to all same value `sigmasq_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#'
#' @return named list of prior parameters
#' @export
set_truncnorm_prior_parameters <- function(
    dims,
    mu_p = sqrt(100/dims$N),
    Mu_p = matrix(mu_p, nrow = dims$K, ncol = dims$N),
    sigmasq_p = mu_p, #mu_p/10,
    Sigmasq_p = matrix(sigmasq_p, nrow = dims$K, ncol = dims$N),
    mu_e = sqrt(100/dims$N), #mean(colSums(M))/(N*100),
    Mu_e = matrix(mu_e, nrow = dims$N, ncol = dims$G),
    sigmasq_e = mu_e, #mu_e/10,
    Sigmasq_e = matrix(sigmasq_e, nrow = dims$N, ncol = dims$G),
    alpha = 0.1,
    Alpha = rep(alpha, dims$K),
    beta = 2,
    Beta = rep(beta, dims$K),
    a = 0.8,
    b = 0.8
) {
    list(
        Mu_p = Mu_p,
        Sigmasq_p = Sigmasq_p,
        Mu_e = Mu_e,
        Sigmasq_e = Sigmasq_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    )
}

#' Set prior parameters for Exponential prior
#'
#' @param dims named list of dimensions N, K, G
#' @param lambda_p see `Lambda_p`
#' @param Lambda_p rate parameter for the exponential prior on `P`, matrix
#' size K x N. Defaults to all same value `lambda_p`
#' @param lambda_e see `Lambda_e`
#' @param Lambda_e rate parameter for the exponential prior on `E`, matrix
#' size N x G. Defaults to all same value `lambda_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#'
#' @return named list of prior parameters
#' @export
set_exponential_prior_parameters <- function(
    dims,
    lambda_p = sqrt(dims$N/100),
    Lambda_p = matrix(lambda_p, nrow = dims$K, ncol = dims$N),
    lambda_e = sqrt(dims$N/100),
    Lambda_e = matrix(lambda_e, nrow = dims$N, ncol = dims$G),
    alpha = 0.1,
    Alpha = rep(alpha, dims$K),
    beta = 2,
    Beta = rep(beta, dims$K),
    a = 0.8,
    b = 0.8
) {
    list(
        Lambda_p = Lambda_p,
        Lambda_e = Lambda_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    )
}

#' Set prior parameters for Gamma prior
#'
#' @param dims named list of dimensions N, K, G
#' @param alpha_p see `Alpha_p`
#' @param Alpha_p shape parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `alpha_p`
#' @param beta_p see `Beta_p`
#' @param Beta_p rate parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `beta_p`
#' @param alpha_e see `Alpha_e`
#' @param Alpha_e shape parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `alpha_e`
#' @param beta_e see `Beta_e`
#' @param Beta_e rate parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `beta_e`
#' #'
#' @return named list of prior parameters
#' @export
set_gamma_prior_parameters <- function(
    dims,
    alpha_p = 10,
    Alpha_p = matrix(alpha_p, nrow = dims$K, ncol = dims$N),
    beta_p = sqrt(dims$N),
    Beta_p = matrix(beta_p, nrow = dims$K, ncol = dims$N),
    alpha_e = 10,
    Alpha_e = matrix(alpha_e, nrow = dims$N, ncol = dims$G),
    beta_e = sqrt(dims$N),
    Beta_e = matrix(beta_e, nrow = dims$N, ncol = dims$G),
    a = 0.8,
    b = 0.8
) {
    list(
        Alpha_p = Alpha_p,
        Beta_p = Beta_p,
        Alpha_e = Alpha_e,
        Beta_e = Beta_e,
        a = a,
        b = b
    )
}

#' initialize Theta
#'
#' @param likelihood string, one of c('poisson','normal')
#' @param prior string, one of c('exponential','truncnormal','gamma')
#' @param learn_A boolean, whether A will be learned or will be fixed
#' @param dims named list of dimensions
#' @param inits list of initial values, optional
#' @param prior_parameters list of named prior parameters, optional
#'
#' @return named list of initialized unknowns
#' @noRd
initialize_Theta <- function(
    likelihood, prior,
    learn_A, dims,
    inits = NULL,
    prior_parameters = NULL
) {
    # prior parameters
    if (is.null(prior_parameters)) {
        if (prior == 'truncnormal') {
            Theta = set_truncnorm_prior_parameters(dims)
        } else if (prior == 'exponential') {
            Theta = set_exponential_prior_parameters(dims)
        } else if (prior == 'gamma') {
            Theta = set_gamma_prior_parameters(dims)
        }
    } else {
        Theta = prior_parameters
    }

    # signatures P
    if (is.null(inits$P)) {
        Theta$P <- matrix(nrow = dims$K, ncol = dims$N)
        for (k in 1:dims$K) {
            for (n in 1:dims$N) {
                if (prior == 'truncnormal') {
                    Theta$P[k,n] <- truncnorm::rtruncnorm(
                        1, a = 0, b = Inf,
                        mean = Theta$Mu_p[k,n],
                        sd = sqrt(Theta$Sigmasq_p[k,n])
                    )
                } else if (prior == 'exponential') {
                    Theta$P[k,n] <- stats::rexp(1, Theta$Lambda_p[k,n])
                } else if (prior == 'gamma') {
                    Theta$P[k,n] <- stats::rgamma(1, Theta$Alpha_p[k,n], Theta$Beta_p[k,n])
                }
            }
        }
    } else {
        Theta$P <- inits$P
    }

    # exposures E
    if (is.null(inits$E)) {
        Theta$E <- matrix(nrow = dims$N, ncol = dims$G)
        for (n in 1:dims$N) {
            for (g in 1:dims$G) {
                if (prior == 'truncnormal') {
                    Theta$E[n,g] <- truncnorm::rtruncnorm(
                        1, a = 0, b = Inf,
                        mean = Theta$Mu_e[n,g],
                        sd = sqrt(Theta$Sigmasq_e[n,g])
                    )
                } else if (prior == 'exponential') {
                    Theta$E[n,g] <- stats::rexp(1, Theta$Lambda_e[n,g])
                } else if (prior == 'gamma') {
                    Theta$E[n,g] <- stats::rgamma(1, Theta$Alpha_e[n,g], Theta$Beta_e[n,g])
                }
            }
        }
    } else {
        Theta$E <- inits$E
    }

    # variance sigmasq
    if (likelihood == 'normal') {
        if (is.null(inits$sigmasq)) {
            Theta$sigmasq <- sapply(1:dims$K, function(k) {
                invgamma::rinvgamma(n = 1, shape = Theta$Alpha[k], scale = Theta$Beta[k])
            })
        } else {
            Theta$sigmasq <- inits$sigmasq
        }
    } else if (likelihood == 'poisson') {
        Theta$Z <- array(dim = c(dims$K, dims$N, dims$G))
        for (k in 1:dims$K) {
            for (g in 1:dims$G) {
                Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
            }
        }
    }

    # signature assignment A
    if (learn_A) {
        Theta$q <- matrix(
            rbeta(dims$S * dims$N, Theta$a, Theta$b),
            nrow = dims$S, ncol = dims$N
        )
        Theta$A <- matrix(
            as.numeric(runif(dims$S * dims$N) < c(Theta$q)),
            nrow = dims$S, ncol = dims$N
        )
    } else if (!is.null(inits$A)) {
        Theta$A <- inits$A
        Theta$q <- inits$A
    } else {
        Theta$A <- matrix(1, nrow = dims$S, ncol = dims$N)
        Theta$q <- Theta$A
    }

    return(Theta)
}

#' Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N number of signatures
#' @param P (optional) initial signatures matrix, K x N
#' @param E (optional) initial exposure matrix, N x G
#' @param sigmasq (optinal) initial variance vector, length K
#' @param prior string, one of c('truncnormal','exponential')
#' @param niters how many iterations the Gibbs sampler is run
#' @param burn_in the first `burn_in` iterations will be ignored when computing
#' MAP estimate
#' @param logevery the log, save, and plot files will be updated every
#' `logevery` iterations
#' @param file file name without extension of log, save, and plot files
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#' @param true_P (optional) true signatures matrix P to compare to in a heatmap
#'
#' @return list
#' @export
bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        inits = NULL,
        likelihood = "normal",
        prior = "truncnormal",
        prior_parameters = NULL,
        logevery = 100,
        file = paste0('nmf_', likelihood, '_', prior),
        overwrite = FALSE,
        true_P = NULL,
        niters = NULL,
        burn_in = NULL
) {

    if (likelihood == 'normal' & prior == 'exponential' & !is.null(max_N)) {
        warning(paste0(
            "Cannot learn N for Normal-Exponential model,",
            "setting N = ", max_N
        ))
        N = max_N
        max_N = NULL
    }

    if (likelihood == 'poisson' & prior == 'exponential' & !is.null(max_N)) {
        warning(paste0(
            "Cannot learn N for Poisson-Exponential model,",
            "setting N = ", max_N
        ))
        N = max_N
        max_N = NULL
    }

    # check N/max_N combination is valid
    if (is.null(N) & is.null(max_N)) {
        stop("Either `N` or `max_N` must be provided.")
    } else if (!is.null(N) & !is.null(max_N)) {
        message("Both `N` and `max_N` provided, using `N`.")
    } else if (is.null(N)) {
        N = max_N
    }

    # check prior and likelihood are valid
    if (!(likelihood %in% c('normal', 'poisson'))) {
        stop("likelihood must be one of c('normal')")
    } else if (likelihood == 'normal') {
        if (!(prior %in% c('truncnormal','exponential'))) {
            stop("prior must be one of c('truncnormal','exponential') with `likelihood = 'normal'`")
        }
    } else if (likelihood == 'poisson') {
        if (!(prior %in% c('gamma','exponential'))) {
            stop("prior must be one of c('gamma','exponential') with `likelihood = 'poisson'`")
        }
    }

    # check whether A is given or learned
    learn_A <- !is.null(max_N) & is.null(inits$A)

    # number of iterations depends on whether we're learning A and likelihood
    if (is.null(niters)) {
        if (!learn_A) {
            niters = 1500
            burn_in = 1000
        } else if (likelihood == 'normal') {
            niters_maxN = 5000
            burn_in_maxN = 3500
        } else if (likleihood == 'poisson') {
            niters_maxN = 10000
            burn_in_maxN = 7500
        }
    }

    # set up tempering schedule
    if (learn_A) {
        gamma_sched <- get_gamma_sched(len = niters)
    } else {
        gamma_sched <- rep(1, niters)
    }

    if (likelihood == 'poisson') {
        logfac = vector(length = max(M))
        logfac[1] = 0
        for (i in 2:length(logfac)) {
            logfac[i] = log(i) + logfac[i-1]
        }
    } else {
        logfac = NULL
    }

    # set up dimensions
    dims = list(
        K = dim(M)[1],
        G = dim(M)[2],
        N = N,
        S = 1
    )

    # check burn_in is valid
    if (burn_in > niters) {
        message(paste0(
            "Burn in ", burn_in, " is greater than niters ",
            niters, ", setting burn_in = 0"
        ))
        burn_in = 0
    }

    # set up file names
    savefile = paste0(file, '.res')
    logfile = paste0(file, '.log')
    plotfile = paste0(file, '.pdf')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        savefile = paste0(file, '_', tail, '.res')
        logfile = paste0(file, '_', tail, '.log')
        plotfile = paste0(file, '_', tail, '.pdf')
    }

    # start logging
    sink(file = logfile)
    print(Sys.time())

    # set up Theta
    Theta <- initialize_Theta(
        likelihood, prior,
        learn_A, dims,
        inits = inits,
        prior_parameters = prior_parameters
    )

    # set up logs
    RMSE <- c()
    KL <- c()
    loglik <- c()

    P.log <- list()
    E.log <- list()
    sigmasq.log <- list()
    A.log <- list()
    q.log <- list()

    # Gibbs sampler: sample niters times
    for (iter in 1:niters) {

        # update P
        for (n in 1:dims$N) {
            Theta$P[, n] <- sample_Pn(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
        }

        # update E
        for (n in 1:dims$N) {
            Theta$E[n, ] <- sample_En(n, M, Theta, dims, likelihood = likelihood, prior = prior, gamma = gamma_sched[iter])
        }

        # if Normal likelihood, update sigmasq
        if (likelihood == 'normal') {
            Theta$sigmasq <- sample_sigmasq_normal(M, Theta, dims, gamma = gamma_sched[iter])
        }

        # if Poisson likelihood, update latent counts Z
        if (likelihood == 'poisson') {
            for (k in 1:dims$K) {
                for (g in 1:dims$G) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
                }
            }
        }

        # update A and q if learn_A = TRUE
        if (learn_A) {
            for (n in 1:dims$N) {
                Theta$A[1, n] <- sample_An(n, M, Theta, dims, logfac, likelihood = likelihood, gamma = gamma_sched[iter])
                Theta$q[1, n] <- sample_qn(n, Theta, gamma = gamma_sched[iter])
            }
        }

        # record in logs
        Mhat <- Theta$P %*% diag(Theta$A[1, ]) %*% Theta$E
        RMSE <- c(RMSE, get_RMSE(M, Mhat))
        KL <- c(KL, get_KLDiv(M, Mhat))
        if (likelihood == 'normal') {
            loglik <- c(loglik, get_loglik_normal(M, Theta, dims))
        } else if (likelihood == 'poisson') {
            loglik <- c(loglik, get_loglik_poisson(M, Theta, dims, logfac))
        }

        P.log[[iter]] <- Theta$P
        E.log[[iter]] <- Theta$E
        sigmasq.log[[iter]] <- Theta$sigmasq
        A.log[[iter]] <- Theta$A
        q.log[[iter]] <- Theta$q

        # periodically save every logevery iters and at the end
        if (iter %% logevery == 0 | iter == niters) {
            cat(paste(iter, "/", niters, "\n"))

            # plot metrics
            grDevices::pdf(plotfile)
            graphics::par(mfrow = c(3,1))
            plot(RMSE)
            if (sum(!is.na(KL)) > 0) {
                if (sum(KL != -Inf, na.rm = TRUE) > 0) {
                    plot(KL)
                }
            }
            if (sum(!is.na(loglik)) > 0){
                if (sum(loglik != -Inf, na.rm = TRUE) > 0) {
                    plot(loglik)
                }
            }
            grDevices::dev.off()

            # only use burn_in if iter > burn_in
            if (iter > burn_in) {
                keep = burn_in:iter
            } else {
                keep = 1:iter
            }

            # get MAP of A matrix (fine to do even if learn_A = FALSE)
            A.map = get_mode(A.log[keep])
            map.idx = keep[A.map$idx]
            if (learn_A) {
                print(A.map$top_counts)
            }

            # only keep signatures present in MAP A matrix
            keep_sigs = as.logical(A.map$matrix[1, ])

            # save results
            res <- list(
                M = M,
                true_P = true_P,
                logs = list(
                    P = P.log,
                    E = E.log,
                    sigmasq = sigmasq.log,
                    q = q.log,
                    A = A.log
                ),
                MAP = list(
                    A = A.map$matrix,
                    q = get_mean(q.log[map.idx]),
                    P = get_mean(P.log[map.idx])[, keep_sigs],
                    E = get_mean(E.log[map.idx])[keep_sigs, ],
                    sigmasq = get_mean(sigmasq.log[map.idx])
                ),
                metrics = list(
                    loglik = loglik,
                    RMSE = RMSE,
                    KL = KL
                ),
                burn_in = burn_in,
                niters = niters,
                final_Theta = Theta,
                dims = dims
            )
            save(res, file = savefile)
        }
    }

    # similarity matrix with true_P, if provided
    # only at end because it slows things down
    if (!is.null(true_P) & dims$N > 1) {
        sim_mat <- pairwise_sim(
            res$MAP$P, true_P,
            name1 = "estimated", name2 = "true",
            which = "cols"
        )
        heatmap <- get_heatmap(res$MAP$P, true_P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
        save(res, file = savefile)
    }

    # end log
    sink()
    return(res)
}
