#' Compute log target pdf for updating En with Poisson-Truncated Normal model
#'
#' @param M mutational catalog matrix, K x G
#' @param En vector length K, value of En to evaluate
#' @param n integer, signature index
#' @param Theta list of parameters
#'
#' @return scalar
#' @noRd
log_target_En_poisson_truncnorm <- function(M, En, n, Theta) {
    Theta$E[n,] <- En
    Mhat <- get_Mhat(Theta)
    log_prior <- log(truncnorm::dtruncnorm(
        En,
        mean = Theta$Mu_e[n,],
        sd = sqrt(Theta$Sigmasq_e[n,]),
        a = 0,
        b = Inf
    ))
    log_likelihood <- colSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' Compute log target pdf for updating En with Poisson-Exponential model
#'
#' @param M mutational catalog matrix, K x G
#' @param En vector length K, value of En to evaluate
#' @param n integer, signature index
#' @param Theta list of parameters
#'
#' @return scalar
#' @noRd
log_target_En_poisson_exp <- function(M, En, n, Theta) {
    Theta$E[n,] <- En
    Mhat <- get_Mhat(Theta)
    log_prior <- log(dexp(
        En, Theta$Lambda_e[n,]
    ))
    log_likelihood <- colSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' Get mu and sigmasq for E[n,] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_exponential <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n], # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_1 <- mu_num_term_1 / Theta$sigmasq
    mu_num_term_2 <- Theta$Lambda_e[n, ] # length G
    denom <- sum(gamma * Theta$A[1,n] * Theta$P[, n] ** 2) / Theta$sigmasq

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
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal_truncnormal <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n], # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G
    mu_num_term_1 <- mu_num_term_1 / Theta$sigmasq
    mu_num_term_2 <- Theta$Mu_e[n, ] / Theta$Sigmasq_e[n,] # length G
    denom <- (1/Theta$Sigmasq_e[n,]) +
        gamma * sum(Theta$A[1,n] * Theta$P[, n] ** 2) / Theta$sigmasq

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
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_En_normal <- function(n, M, Theta, dims, prior, gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_truncnormal(
            n, M, Theta, dims, gamma = gamma
        )
    } else if (prior == 'exponential') {
        if (Theta$A[1,n] == 0) {
            # sample from prior, doesn't collapse like truncnormal
            sampled <- stats::rexp(dims$G, Theta$Lambda_e[n,])
            return(sampled)
        }
        mu_sigmasq_E <- get_mu_sigmasq_En_normal_exponential(
            n, M, Theta, dims, gamma = gamma
        )
    }
    mu_E = mu_sigmasq_E$mu
    sigmasq_E = mu_sigmasq_E$sigmasq

    # sample from truncated normal
    sampled <- truncnorm::rtruncnorm(1, mean = mu_E, sd = sqrt(sigmasq_E), a = 0, b = Inf)
    return(sampled)
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
sample_En_poisson <- function(n, M, Theta, dims, prior, gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(
                1,
                Theta$Alpha_e[n,g] + gamma*sum(Theta$Z[,n,g]),
                Theta$Beta_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n])
            )
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$G, function(g) {
            rgamma(
                1,
                1 + gamma * sum(Theta$Z[,n,g]),
                Theta$Lambda_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n])
            )
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
#' @param gamma double, tempering parameter. DO NOT CHANGE, KEEP gamma = 1
#'
#' @return vector of length G
#' @noRd
sample_En <- function(n, M, Theta, dims, likelihood, prior, fast, gamma = 1) {
    if (likelihood == 'normal' | (likelihood == 'poisson' & fast)) {
        proposal <- sample_En_normal(n, M, Theta, dims, prior = prior, gamma = gamma)
        if (fast) {
            if (prior == 'truncnormal') {
                acceptance <- exp(
                    log_target_En_poisson_truncnorm(M, proposal, n, Theta) -
                    log_target_En_poisson_truncnorm(M, Theta$E[n,], n, Theta)
                )
            } else {
                acceptance <- exp(
                    log_target_En_poisson_exp(M, proposal, n, Theta) -
                    log_target_En_poisson_exp(M, Theta$E[n,], n, Theta)
                )
            }
            # acceptance prob is NaN if 0/0, so give it a 50-50 chance
            acceptance[is.na(acceptance)] <- 0.5
            accepted <- runif(dims$G) < acceptance
            sampled <- Theta$E[n,]
            sampled[accepted] <- proposal[accepted]
        } else {
            sampled <- proposal
            acceptance <- 1
        }
    } else {
        # likelihood == 'poisson' & !fast
        sampled <- sample_En_poisson(n, M, Theta, dims, prior, gamma = gamma)
    }
    return(list(
        sampled = sampled,
        acceptance = acceptance
    ))
}
