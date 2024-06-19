#' Compute log target pdf for updating Pn with Poisson-Truncated Normal model
#'
#' @param M mutational catalog matrix, K x G
#' @param Pn vector length K, value of Pn to evaluate
#' @param n integer, signature index
#' @param Theta list of parameters
#'
#' @return scalar
#' @noRd
log_target_Pn_poisson_truncnorm <- function(M, Pn, n, Theta) {
    Theta$P[,n] <- Pn
    Mhat <- get_Mhat(Theta)
    log_prior <- log(truncnorm::dtruncnorm(
        Pn,
        mean = Theta$Mu_p[,n],
        sd = sqrt(Theta$Sigmasq_p[,n]),
        a = 0,
        b = Inf
    ))
    log_likelihood <- rowSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' #' Compute log proposal pdf for updating Pn with Normal-Truncated Normal model
#' #'
#' #' @param M mutational catalog matrix, K x G
#' #' @param Pn vector length K, value of Pn to evaluate
#' #' @param n integer, signature index
#' #' @param Theta list of parameters
#' #'
#' #' @return scalar
#' #' @noRd
#' log_proposal_Pn_normal_truncnorm <- function(M, Pn, n, Theta) {
#'     Theta$P[,n] <- Pn
#'     Mhat <- get_Mhat(Theta)
#'     log_prior <- log(truncnorm::dtruncnorm(
#'         Pn,
#'         mean = Theta$Mu_p[,n],
#'         sd = sqrt(Theta$Sigmasq_p[,n]),
#'         a = 0,
#'         b = Inf
#'     ))
#'     sigmasq_mat <- matrix(
#'         rep(Theta$sigmasq, nrow(M)),
#'         nrow = nrow(M),
#'         byrow = TRUE
#'     )
#'     log_likelihood <- rowSums(dnorm(M, mean = Mhat, sd = sqrt(sigmasq_mat), log = TRUE))
#'     log_proposal <- log_prior + log_likelihood
#'     return(log_proposal)
#' }

#' Compute log target pdf for updating Pn with Poisson-Exponential model
#'
#' @param M mutational catalog matrix, K x G
#' @param Pn vector length K, value of Pn to evaluate
#' @param n integer, signature index
#' @param Theta list of parameters
#'
#' @return scalar
#' @noRd
log_target_Pn_poisson_exp <- function(M, Pn, n, Theta) {
    Theta$P[,n] <- Pn
    Mhat <- get_Mhat(Theta)
    log_prior <- log(dexp(
        Pn, Theta$Lambda_p[,n]
    ))
    log_likelihood <- rowSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' #' Compute log proposal pdf for updating Pn with Normal-Exponential model
#' #'
#' #' @param M mutational catalog matrix, K x G
#' #' @param Pn vector length K, value of Pn to evaluate
#' #' @param n integer, signature index
#' #' @param Theta list of parameters
#' #'
#' #' @return scalar
#' #' @noRd
#' log_proposal_Pn_normal_exp <- function(M, Pn, n, Theta) {
#'     Theta$P[,n] <- Pn
#'     Mhat <- get_Mhat(Theta)
#'     log_prior <- log(dexp(
#'         Pn, Theta$Lambda_p[,n]
#'     ))
#'     sigmasq_mat <- matrix(
#'         rep(Theta$sigmasq, nrow(M)),
#'         nrow = nrow(M),
#'         byrow = TRUE
#'     )
#'     log_likelihood <- rowSums(dnorm(M, mean = Mhat, sd = sqrt(sigmasq_mat), log = TRUE))
#'     log_proposal <- log_prior + log_likelihood
#'     return(log_proposal)
#' }

#' get mu and sigmasq for P[,n] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_exponential <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * Theta$A[1,n] * sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$E[n, ] / Theta$sigmasq, # length G
        "*"
    ) %>% # dim KxG
        rowSums() # length K
    mu_num_term_2 <- Theta$Lambda_p[, n] # length K
    denom <- gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2 / Theta$sigmasq)

    mu_P <- (mu_num_term_1 - mu_num_term_2) / denom # length K
    sigmasq_P <- 1 / denom # length K
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
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_truncnormal <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    # compute mean
    mu_num_term_1 <- gamma * (sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$A[1,n] * Theta$E[n, ] / Theta$sigmasq, # length G
        "*"
    ) %>% # dim KxG
        rowSums()) # length K
    mu_num_term_2 <- Theta$Mu_p[, n] / Theta$Sigmasq_p[,n] # length K
    denom <- (1/Theta$Sigmasq_p[,n]) +
        gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2 / Theta$sigmasq)


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
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_normal <- function(n, M, Theta, dims, prior, gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_truncnormal(
            n, M, Theta, dims, gamma = gamma
        )
    } else if (prior == 'exponential') {
        if (Theta$A[1,n] == 0) {
            # sample from prior, doesn't collapse like truncnormal
            sampled <- stats::rexp(dims$K, Theta$Lambda_p[,n])
            return(sampled)
        }
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_exponential(
            n, M, Theta, dims, gamma = gamma
        )
    }

    mu_P = mu_sigmasq_P$mu
    sigmasq_P = mu_sigmasq_P$sigmasq

    # sample from truncated normal
    sampled <- truncnorm::rtruncnorm(
        1, mean = mu_P, sd = sqrt(sigmasq_P), a = 0, b = Inf
    )
    return(sampled)
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
sample_Pn_poisson <- function(n, M, Theta, dims, prior, gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(
                1,
                Theta$Alpha_p[k,n] + gamma * sum(Theta$Z[k,n,]),
                Theta$Beta_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,])
            )
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(
                1,
                1 + gamma * sum(Theta$Z[k,n,]),
                Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,])
            )
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
#' @param gamma double, tempering parameter. DO NOT CHANGE, KEEP gamma = 1
#'
#' @return vector length K
#' @noRd
sample_Pn <- function(n, M, Theta, dims, likelihood, prior, fast, gamma = 1) {
    if (likelihood == 'normal' | (likelihood == 'poisson' & fast)) {
        sampled <- sample_Pn_normal(n, M, Theta, dims, prior, gamma)
        # if (fast) {
        #     if (prior == "truncnormal") {
        #         acceptance <- exp(
        #             log_target_Pn_poisson_truncnorm(M, proposal, n, Theta) +
        #             log_proposal_Pn_normal_truncnorm(M, Theta$P[,n], n, Theta) -
        #             log_target_Pn_poisson_truncnorm(M, Theta$P[,n], n, Theta) -
        #             log_proposal_Pn_normal_truncnorm(M, proposal, n, Theta)
        #         )
        #     } else {
        #         acceptance <- exp(
        #             log_target_Pn_poisson_exp(M, proposal, n, Theta) +
        #             log_proposal_Pn_normal_exp(M, Theta$P[,n], n, Theta) -
        #             log_target_Pn_poisson_exp(M, Theta$P[,n], n, Theta) -
        #             log_proposal_Pn_normal_exp(M, proposal, n, Theta)
        #         )
        #     }
        #     # acceptance prob is NaN if 0/0, so give it a 50-50 chance
        #     acceptance[is.na(acceptance)] <- 0.5
        #     accepted <- runif(dims$K) < acceptance
        #     sampled <- Theta$P[,n]
        #     sampled[accepted] <- proposal[accepted]
        # } else {
        #     sampled <- proposal
        #     acceptance <- 1
        # }
    } else {
        # likelihood == 'poisson' & !fast
        sampled <- sample_Pn_poisson(n, M, Theta, dims, prior, gamma)
    }
    acceptance <- 1
    return(list(
        sampled = sampled,
        acceptance = acceptance
    ))
}
