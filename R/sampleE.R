#' Compute log target pdf for updating En with Poisson model
#'
#' @param M mutational catalog matrix, K x G
#' @param En vector length K, value of En to evaluate
#' @param n integer, signature index
#' @param Theta list, named list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#'
#' @return scalar
#' @noRd
log_target_En <- function(M, En, n, Theta, prior) {
    Theta$E[n,] <- En
    Mhat <- get_Mhat(Theta)
    if (prior == 'truncnormal') {
        log_prior <- log(truncnorm::dtruncnorm(
            En, mean = Theta$Mu_e[n,], sd = sqrt(Theta$Sigmasq_e[n,]),
            a = 0, b = Inf
        ))
    } else if (prior == 'exponential') {
        log_prior <- dexp(
            En, Theta$Lambda_e[n,], log = TRUE
        )
    }
    log_likelihood <- colSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' Get mu and sigmasq for E[n,] in Normal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_En_normal <- function(n, M, Theta, prior, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, n)

    # broadcast Pn to residuals M - Mhat
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        1, # multiply each column by P[,n]
        Theta$P[, n], # length K
        "*"
    ) %>% # dim KxG
        colSums() # length G

    mu_num_term_1 <- mu_num_term_1 / Theta$sigmasq
    denom <- sum(gamma * Theta$A[1,n] * Theta$P[, n] ** 2) / Theta$sigmasq

    if (prior == 'exponential') {
        mu_num_term_2 <- Theta$Lambda_e[n, ] # length G
        mu_E <- (mu_num_term_1 - mu_num_term_2) / denom # length G
        sigmasq_E <- 1 / denom
    } else if (prior == 'truncnormal') {
        mu_num_term_2 <- Theta$Mu_e[n, ] / Theta$Sigmasq_e[n,] # length G
        denom <- denom + (1/Theta$Sigmasq_e[n,])
        mu_E <- (mu_num_term_1 + mu_num_term_2) / denom # length G
        sigmasq_E <- 1 / denom
    }

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

    # Normal-exponential doesn't collapse like normal-truncated normal
    # sample from prior when Ann = 0 or gamma = 0
    if (prior == 'exponential' & (Theta$A[1,n] == 0 | gamma == 0)) {
        sampled <- stats::rexp(dims$G, Theta$Lambda_e[n,])
        return(sampled)
    }

    # Otherwise, compute mean and sd for sample from full conditional
    mu_sigmasq_E <- get_mu_sigmasq_En_normal(n, M, Theta, prior, gamma = gamma)

    # Sample from truncated normal
    sampled <- truncnorm::rtruncnorm(
        1, mean = mu_sigmasq_E$mu, sd = sqrt(mu_sigmasq_E$sigmasq),
        a = 0, b = Inf
    )
    return(sampled)
}

#' Sample E[n,] for Poisson likelihood
#'
#' @param n signature index
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
        shape <- sapply(1:dims$G, function(g) {
            Theta$Alpha_e[n,g] + gamma*sum(Theta$Z[,n,g])
        })
        rate <- sapply(1:dims$G, function(g) {
            Theta$Beta_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n])
        })
    } else if (prior == 'exponential') {
        shape <- sapply(1:dims$G, function(g) {
            1 + gamma * sum(Theta$Z[,n,g])
        })
        rate <- sapply(1:dims$G, function(g) {
            Theta$Lambda_e[n,g] + gamma * Theta$A[1,n] * sum(Theta$P[,n])
        })
    }
    sampled <- sapply(1:dims$G, function(g) { rgamma(1, shape[g], rate[g]) })
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
        proposal <- sample_En_normal(n, M, Theta, dims, prior, gamma)
        if (likelihood == 'poisson' & fast) {
            acceptance <- exp(
                log_target_En(M, proposal, n, Theta, prior) -
                log_target_En(M, Theta$E[n,], n, Theta, prior)
            )

            acceptance[is.na(acceptance)] <- 0.5 # NaN if 0/0, give 50-50 chance
            acceptance[acceptance > 1] <- 1 # cap acceptance probability at 1

            # accept with probability acceptance
            accepted <- runif(dims$G) < acceptance
            sampled <- Theta$E[n,]
            sampled[accepted] <- proposal[accepted]
        } else {
            # likelihood == 'normal'
            sampled <- proposal
            acceptance <- 1
        }
    } else {
        # likelihood == 'poisson' & !fast
        sampled <- sample_En_poisson(n, M, Theta, dims, prior, gamma = gamma)
        acceptance <- 1
    }
    return(list(
        sampled = sampled,
        acceptance = acceptance
    ))
}
